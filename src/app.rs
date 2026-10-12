use crate::theme::{self, Palette, stroke};
use dna_analyzer::alignment::{self, AlignmentResult};
use dna_analyzer::restriction::{self, RestrictionHit};
use dna_analyzer::{
    MoleculeKind, SequenceAnalysis, analyze_input, find_overlapping, format_origin_export,
    six_frame_translation,
};
use eframe::egui::{
    self, Align, Color32, FontId, Label, Layout, Margin, RichText, ScrollArea, TextEdit,
    TextFormat, TextStyle, text::LayoutJob,
};
use std::fs;
use std::io::Write;
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::time::{Duration, Instant};

const BUTTON_HEIGHT: f32 = 28.0;
const INPUT_HEIGHT: f32 = 145.0;
const OUTPUT_HEIGHT: f32 = 92.0;

struct RestrictionTask {
    receiver: Receiver<Vec<RestrictionHit>>,
    started: Instant,
}

struct AlignmentTask {
    receiver: Receiver<Result<AlignmentResult, String>>,
    started: Instant,
}

/// A restriction hit with its display text prepared once, not every frame.
struct RestrictionRow {
    name: String,
    site: String,
    site_label: String,
    positions: String,
}

impl From<RestrictionHit> for RestrictionRow {
    fn from(hit: RestrictionHit) -> Self {
        let positions = hit
            .positions
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        Self {
            site_label: format!(" ({}): ", hit.site),
            name: hit.name,
            site: hit.site,
            positions,
        }
    }
}

impl RestrictionRow {
    fn matches(&self, keywords: &[String]) -> bool {
        keywords.is_empty()
            || keywords.iter().any(|keyword| {
                contains_ignore_ascii_case(&self.name, keyword)
                    || contains_ignore_ascii_case(&self.site, keyword)
            })
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AppMode {
    Main,
    Search,
    Alignment,
}

#[derive(Clone, Copy)]
enum MainAction {
    Clear,
    Translate,
    Restriction,
    Blast,
    Alignment,
    Save,
}

const MAIN_ACTIONS: [(&str, MainAction); 6] = [
    ("Clear", MainAction::Clear),
    ("Translate 6-Frames", MainAction::Translate),
    ("Show Restriction Enzymes", MainAction::Restriction),
    ("NCBI BLAST", MainAction::Blast),
    ("Multiple Alignment", MainAction::Alignment),
    ("Save as Fasta", MainAction::Save),
];

pub struct DnaAnalyzerApp {
    mode: AppMode,
    input: String,
    fasta_name: String,
    analysis: SequenceAnalysis,
    translations: [String; 6],
    translation_visible: bool,
    status: String,
    status_is_error: bool,

    restriction_visible: bool,
    restriction_filter: String,
    restriction_rows: Vec<RestrictionRow>,
    restriction_task: Option<RestrictionTask>,
    /// Whether `restriction_rows` belongs to the current sequence.
    restriction_fresh: bool,

    search_query: String,
    search_query_len: usize,
    search_matches: Vec<usize>,

    alignment_inputs: Vec<String>,
    alignment_result: Option<AlignmentResult>,
    alignment_task: Option<AlignmentTask>,
    child_processes: Vec<Child>,
}

impl DnaAnalyzerApp {
    pub fn new(
        creation_context: &eframe::CreationContext<'_>,
        mode: AppMode,
        initial_input: Option<String>,
    ) -> Self {
        theme::install(&creation_context.egui_ctx);

        let mut app = Self {
            mode,
            input: String::new(),
            fasta_name: "sequence_1".to_owned(),
            analysis: SequenceAnalysis::default(),
            translations: std::array::from_fn(|_| String::new()),
            translation_visible: false,
            status: String::new(),
            status_is_error: false,
            restriction_visible: false,
            restriction_filter: String::new(),
            restriction_rows: Vec::new(),
            restriction_task: None,
            restriction_fresh: false,
            search_query: String::new(),
            search_query_len: 0,
            search_matches: Vec::new(),
            alignment_inputs: vec![String::new(), String::new()],
            alignment_result: None,
            alignment_task: None,
            child_processes: Vec::new(),
        };
        if let Some(input) = initial_input {
            app.input = input;
            app.refresh_analysis();
        }
        if mode != AppMode::Main {
            return app;
        }

        let mut startup_file = None;
        let mut show_translation = false;
        let mut show_restriction = false;
        for argument in std::env::args_os().skip(1) {
            match argument.to_string_lossy().as_ref() {
                "--show-translation" => show_translation = true,
                "--show-restriction" => show_restriction = true,
                _ => {
                    let path = std::path::PathBuf::from(argument);
                    if path.is_file() {
                        startup_file = Some(path);
                    }
                }
            }
        }
        if let Some(path) = startup_file {
            app.load_file(&path);
        }
        app.translation_visible = show_translation;
        if show_restriction && !app.analysis.dna.is_empty() {
            app.restriction_visible = true;
            app.start_restriction_analysis();
        }
        app
    }

    fn refresh_analysis(&mut self) {
        self.analysis = analyze_input(&self.input);
        self.translations = six_frame_translation(&self.analysis.dna);
        self.refresh_search();
        self.restriction_task = None;
        self.restriction_rows.clear();
        self.restriction_fresh = false;
        if self.input.trim().is_empty() {
            self.set_status("");
        } else if self.analysis.dna.is_empty() {
            self.set_error("No valid DNA or RNA bases found.");
        } else {
            self.set_status(format!(
                "Parsed {} {} bases ({}).",
                self.analysis.dna.len(),
                self.analysis.kind.label(),
                self.analysis.input_format.label()
            ));
        }
    }

    fn refresh_search(&mut self) {
        self.search_query_len = analyze_input(&self.search_query).dna.len();
        self.search_matches = find_overlapping(&self.analysis.dna, &self.search_query);
    }

    fn poll_tasks(&mut self, context: &egui::Context) {
        let restriction_update =
            self.restriction_task
                .as_ref()
                .and_then(|task| match task.receiver.try_recv() {
                    Ok(hits) => Some(Ok((hits, task.started.elapsed()))),
                    Err(TryRecvError::Disconnected) => {
                        Some(Err("Restriction analysis worker stopped unexpectedly."))
                    }
                    Err(TryRecvError::Empty) => None,
                });
        if let Some(update) = restriction_update {
            self.restriction_task = None;
            match update {
                Ok((hits, elapsed)) => {
                    self.set_status(format!(
                        "Found {} cutting enzymes in {:.1} ms.",
                        hits.len(),
                        elapsed.as_secs_f64() * 1_000.0
                    ));
                    self.restriction_rows = hits.into_iter().map(RestrictionRow::from).collect();
                    self.restriction_fresh = true;
                }
                Err(message) => self.set_error(message),
            }
        }

        let alignment_update =
            self.alignment_task
                .as_ref()
                .and_then(|task| match task.receiver.try_recv() {
                    Ok(result) => Some((result, task.started.elapsed())),
                    Err(TryRecvError::Disconnected) => Some((
                        Err("Alignment worker stopped unexpectedly.".to_owned()),
                        task.started.elapsed(),
                    )),
                    Err(TryRecvError::Empty) => None,
                });
        if let Some((result, elapsed)) = alignment_update {
            self.alignment_task = None;
            match result {
                Ok(alignment) => {
                    self.set_status(format!(
                        "Aligned {} sequences with {} in {:.2} s.",
                        alignment.sequences.len(),
                        alignment.engine,
                        elapsed.as_secs_f64()
                    ));
                    self.alignment_result = Some(alignment);
                }
                Err(message) => self.set_error(message),
            }
        }

        if self.restriction_task.is_some() || self.alignment_task.is_some() {
            context.request_repaint_after(Duration::from_millis(50));
        }

        self.child_processes
            .retain_mut(|child| matches!(child.try_wait(), Ok(None)));
    }

    fn handle_dropped_files(&mut self, context: &egui::Context) {
        let path = context.input(|input| {
            input
                .raw
                .dropped_files
                .iter()
                .find_map(|file| file.path.clone())
        });
        if let Some(path) = path {
            self.load_file(&path);
        }
    }

    fn set_status(&mut self, message: impl Into<String>) {
        self.status = message.into();
        self.status_is_error = false;
    }

    fn set_error(&mut self, message: impl Into<String>) {
        self.status = message.into();
        self.status_is_error = true;
    }

    fn load_file(&mut self, path: &std::path::Path) {
        match fs::read_to_string(path) {
            Ok(contents) => {
                self.input = contents;
                if let Some(stem) = path.file_stem().and_then(|value| value.to_str()) {
                    self.fasta_name = stem.to_owned();
                }
                self.refresh_analysis();
                if !self.analysis.dna.is_empty() {
                    self.set_status(format!("Loaded {} bases.", self.analysis.dna.len()));
                }
            }
            Err(error) => self.set_error(format!("Could not open {}: {error}", path.display())),
        }
    }

    fn save_origin(&mut self) {
        if self.analysis.dna.is_empty() {
            self.set_error("Enter a valid sequence before saving.");
            return;
        }
        let suggested = format!("{}.fasta", self.fasta_name.trim().replace(' ', "_"));
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name(&suggested)
            .add_filter("Fasta", &["fasta", "fa"])
            .save_file()
        {
            let contents = format_origin_export(&self.fasta_name, &self.analysis.display);
            match fs::write(&path, contents) {
                Ok(()) => {
                    self.set_status(format!(
                        "Saved numbered ORIGIN sequence to {}.",
                        path.display()
                    ));
                }
                Err(error) => self.set_error(format!("Could not save {}: {error}", path.display())),
            }
        }
    }

    fn open_blast(&mut self) {
        if self.analysis.dna.is_empty() {
            self.set_error("Enter sequence first.");
            return;
        }
        let url = format!(
            "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={}",
            self.analysis.dna
        );
        match webbrowser::open(&url) {
            Ok(()) => self.set_status("Opened NCBI BLAST in the default browser."),
            Err(error) => self.set_error(format!("Could not open the browser: {error}")),
        }
    }

    fn start_restriction_analysis(&mut self) {
        if self.analysis.dna.is_empty() {
            self.set_error("Enter sequence first.");
            return;
        }
        let sequence = self.analysis.dna.clone();
        let bases = sequence.len();
        let (sender, receiver) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = sender.send(restriction::analyze_restriction_sites(&sequence));
        });
        self.restriction_task = Some(RestrictionTask {
            receiver,
            started: Instant::now(),
        });
        self.set_status(format!(
            "Scanning {bases} bases against {} enzymes…",
            restriction::enzymes().len()
        ));
    }

    fn start_alignment(&mut self) {
        if self.alignment_task.is_some() {
            return;
        }
        let sequences = self.alignment_inputs.clone();
        let (sender, receiver) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = sender.send(alignment::run_alignment(&sequences));
        });
        self.alignment_task = Some(AlignmentTask {
            receiver,
            started: Instant::now(),
        });
        self.set_status("Running multiple sequence alignment…");
    }

    fn open_search_window(&mut self) {
        if self.analysis.dna.is_empty() {
            self.set_error("Enter sequence first.");
            return;
        }
        let spawned = spawn_window(
            "--search-window",
            "sequence search",
            Some(&self.analysis.display),
        );
        self.track_window("sequence search", spawned);
    }

    fn track_window(&mut self, name: &str, spawned: Result<Child, String>) {
        match spawned {
            Ok(child) => {
                self.child_processes.push(child);
                self.set_status(format!("Opened {name} in a new window."));
            }
            Err(message) => self.set_error(message),
        }
    }

    /// Status or input warnings on the left, sequence statistics on the right.
    fn status_line(&self, ui: &mut egui::Ui, show_stats: bool) {
        let palette = Palette::of(ui);
        let mut warnings = Vec::new();
        match self.analysis.kind {
            MoleculeKind::Rna => warnings.push("RNA detected, using A:U base pairing.".to_owned()),
            MoleculeKind::Mixed => warnings
                .push("Both T and U detected; U is normalized to T for analysis.".to_owned()),
            MoleculeKind::Dna => {}
        }
        if !self.analysis.invalid_letters.is_empty() {
            let invalid: String = self.analysis.invalid_letters.iter().collect();
            warnings.push(format!("Invalid characters removed: {invalid}"));
        }
        let (message, color) = if !warnings.is_empty() {
            (warnings.join("  ·  "), palette.warn)
        } else if self.status_is_error {
            (self.status.clone(), palette.error)
        } else {
            (self.status.clone(), palette.text_weak)
        };

        ui.horizontal(|ui| {
            ui.label(RichText::new(message).color(color));
            if show_stats && !self.analysis.dna.is_empty() {
                ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                    let stats = format!(
                        "{}  ·  {} bp  ·  GC {:.1}%",
                        self.analysis.kind.label(),
                        group_thousands(self.analysis.dna.len()),
                        self.analysis.gc_percent
                    );
                    ui.label(RichText::new(stats).color(palette.text_weak));
                });
            }
        });
    }

    fn main_analysis_ui(&mut self, ui: &mut egui::Ui) {
        self.status_line(ui, true);
        theme::caption(ui, "Enter DNA or RNA sequence (A/T/G/C/U):");
        let response = ScrollArea::vertical()
            .id_salt("main_sequence_input")
            .max_height(INPUT_HEIGHT)
            .min_scrolled_height(INPUT_HEIGHT)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                ui.add_sized(
                    [ui.available_width(), INPUT_HEIGHT],
                    TextEdit::multiline(&mut self.input)
                        .font(TextStyle::Monospace)
                        .hint_text("Paste DNA, RNA, FASTA, or GenBank ORIGIN text here…")
                        .margin(Margin::symmetric(8, 6))
                        .lock_focus(true),
                )
            })
            .inner;
        if response.changed() {
            self.refresh_analysis();
        }

        let has_sequence = !self.analysis.dna.is_empty();
        if theme::wide_button(ui, has_sequence, true, "Search Subsequence").clicked() {
            self.open_search_window();
        }

        let mut action = None;
        ui.columns(MAIN_ACTIONS.len(), |columns| {
            for (column, (label, candidate)) in columns.iter_mut().zip(MAIN_ACTIONS) {
                if theme::wide_button(column, true, false, label).clicked() {
                    action = Some(candidate);
                }
            }
        });
        match action {
            Some(MainAction::Clear) => {
                self.input.clear();
                self.refresh_analysis();
            }
            Some(MainAction::Translate) => self.translation_visible = true,
            Some(MainAction::Restriction) => {
                self.restriction_visible = true;
                self.start_restriction_analysis();
            }
            Some(MainAction::Blast) => self.open_blast(),
            Some(MainAction::Alignment) => {
                let spawned = spawn_window("--alignment-window", "multiple alignment", None);
                self.track_window("multiple alignment", spawned);
            }
            Some(MainAction::Save) => self.save_origin(),
            None => {}
        }

        ui.add_space(6.0);
        let outputs = [
            ("Complement:", "Copy Complement", &self.analysis.complement),
            ("Reverse:", "Copy Reverse", &self.analysis.reverse),
            (
                "Reverse-Complement:",
                "Copy Reverse-Complement",
                &self.analysis.reverse_complement,
            ),
        ];
        let mut copied = None;
        ui.columns(outputs.len(), |columns| {
            for (column, (title, copy_label, sequence)) in columns.iter_mut().zip(outputs) {
                theme::caption(column, title);
                output_box(column, title, sequence);
                if theme::wide_button(column, !sequence.is_empty(), false, copy_label).clicked() {
                    column.ctx().copy_text(sequence.clone());
                    copied = Some(title.trim_end_matches(':'));
                }
            }
        });
        if let Some(title) = copied {
            self.set_status(format!("Copied {} to the clipboard.", title.to_lowercase()));
        }

        if self.translation_visible {
            ui.add_space(6.0);
            for (row, frames) in self.translations.chunks(3).enumerate() {
                ui.columns(frames.len(), |columns| {
                    for (index, (column, protein)) in columns.iter_mut().zip(frames).enumerate() {
                        let title = format!("Frame {}:", row * 3 + index + 1);
                        theme::caption(column, &title);
                        output_box(column, &title, protein);
                    }
                });
            }
        }
    }

    fn restriction_panel_ui(&mut self, ui: &mut egui::Ui) {
        let palette = Palette::of(ui);
        ui.label(RichText::new("Restriction Enzyme Sites (1-based)").strong());
        ui.add(
            TextEdit::singleline(&mut self.restriction_filter)
                .hint_text("Filter enzymes (comma-separated)")
                .margin(Margin::symmetric(8, 5))
                .desired_width(f32::INFINITY),
        );
        if self.restriction_task.is_some() {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(RichText::new("Analyzing…").color(palette.text_weak));
            });
        }
        let needs_refresh = !self.analysis.dna.is_empty()
            && self.restriction_task.is_none()
            && !self.restriction_fresh;
        if needs_refresh {
            ui.label(
                RichText::new("Sequence changed; click Show Restriction Enzymes to refresh.")
                    .color(palette.warn),
            );
        }

        let keywords: Vec<String> = self
            .restriction_filter
            .split(',')
            .map(str::trim)
            .filter(|word| !word.is_empty())
            .map(str::to_owned)
            .collect();
        let font = TextStyle::Monospace.resolve(ui.style());
        let name_format = text_format(&font, palette.accent);
        let site_format = text_format(&font, palette.text_weak);
        let positions_format = text_format(&font, palette.text);
        let mut job = LayoutJob::default();
        for row in self
            .restriction_rows
            .iter()
            .filter(|row| row.matches(&keywords))
        {
            if !job.is_empty() {
                job.append("\n", 0.0, positions_format.clone());
            }
            job.append(&row.name, 0.0, name_format.clone());
            job.append(&row.site_label, 0.0, site_format.clone());
            job.append(&row.positions, 0.0, positions_format.clone());
        }
        let placeholder = if !job.is_empty() || self.restriction_task.is_some() {
            None
        } else if self.analysis.dna.is_empty() {
            Some("Enter a sequence to analyze.")
        } else if needs_refresh {
            Some("Restriction results are waiting to be refreshed.")
        } else if self.restriction_rows.is_empty() {
            Some("No cutting enzymes found.")
        } else {
            Some("No enzymes match the filter.")
        };

        theme::surface_frame(ui).show(ui, |ui| {
            ui.set_min_size(ui.available_size());
            ScrollArea::vertical()
                .id_salt("restriction_results")
                .auto_shrink([false, false])
                .show(ui, |ui| match placeholder {
                    Some(text) => {
                        ui.label(RichText::new(text).color(palette.text_weak));
                    }
                    None => {
                        ui.add(Label::new(job).wrap());
                    }
                });
        });
    }

    fn search_window_ui(&mut self, ui: &mut egui::Ui) {
        let palette = Palette::of(ui);
        let response = ui.add(
            TextEdit::singleline(&mut self.search_query)
                .font(TextStyle::Monospace)
                .hint_text("Search sequence…")
                .margin(Margin::symmetric(8, 6))
                .desired_width(f32::INFINITY),
        );
        if response.changed() {
            self.refresh_search();
        }
        let summary = match (self.search_query_len, self.search_matches.len()) {
            (0, _) => format!(
                "{} bases. Type a subsequence to highlight every occurrence.",
                group_thousands(self.analysis.dna.len())
            ),
            (_, 1) => "1 match".to_owned(),
            (_, count) => format!("{} matches", group_thousands(count)),
        };
        theme::caption(ui, &summary);

        let layout = origin_layout(
            &self.analysis.display,
            &self.search_matches,
            self.search_query_len,
            palette,
            TextStyle::Monospace.resolve(ui.style()),
        );
        theme::surface_frame(ui).show(ui, |ui| {
            ui.set_min_size(ui.available_size());
            ScrollArea::both()
                .auto_shrink([false, false])
                .show(ui, |ui| ui.add(Label::new(layout).extend()));
        });
    }

    fn alignment_window_ui(&mut self, ui: &mut egui::Ui) {
        let palette = Palette::of(ui);
        self.status_line(ui, false);
        ScrollArea::vertical()
            .id_salt("alignment_sequence_inputs")
            .max_height(310.0)
            .show(ui, |ui| {
                for (index, sequence) in self.alignment_inputs.iter_mut().enumerate() {
                    ui.add_sized(
                        [ui.available_width(), 58.0],
                        TextEdit::multiline(sequence)
                            .font(TextStyle::Monospace)
                            .margin(Margin::symmetric(8, 6))
                            .hint_text(format!(
                                "Sequence {} — enter DNA/RNA or protein sequence…",
                                index + 1
                            )),
                    );
                }
            });

        let can_remove = self.alignment_inputs.len() > 2;
        ui.columns(2, |columns| {
            if theme::wide_button(&mut columns[0], true, false, "Add Sequence").clicked() {
                self.alignment_inputs.push(String::new());
            }
            if theme::wide_button(&mut columns[1], can_remove, false, "Remove Sequence").clicked() {
                self.alignment_inputs.pop();
            }
        });

        let running = self.alignment_task.is_some();
        let result_height = (ui.available_height() - BUTTON_HEIGHT - 36.0).max(120.0);
        if let Some(result) = &self.alignment_result {
            theme::caption(ui, &result.engine);
            let layout =
                alignment_layout(result, palette, TextStyle::Monospace.resolve(ui.style()));
            theme::surface_frame(ui).show(ui, |ui| {
                ui.set_min_width(ui.available_width());
                ScrollArea::both()
                    .id_salt("alignment_result")
                    .max_height(result_height)
                    .min_scrolled_height(result_height)
                    .auto_shrink([false, false])
                    .show(ui, |ui| ui.add(Label::new(layout).extend()));
            });
        } else {
            theme::caption(ui, "Alignment");
            theme::surface_frame(ui).show(ui, |ui| {
                ui.set_min_size(egui::vec2(ui.available_width(), result_height));
                ui.centered_and_justified(|ui| {
                    if running {
                        ui.spinner();
                    } else {
                        ui.label(
                            RichText::new("Aligned sequences will appear here.")
                                .color(palette.text_weak),
                        );
                    }
                });
            });
        }

        let label = if running { "Aligning…" } else { "Align" };
        if theme::wide_button(ui, !running, true, label).clicked() {
            self.start_alignment();
        }
    }
}

impl Drop for DnaAnalyzerApp {
    fn drop(&mut self) {
        for child in &mut self.child_processes {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

impl eframe::App for DnaAnalyzerApp {
    fn update(&mut self, context: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_tasks(context);
        let style = context.style();
        let central = egui::CentralPanel::default().frame(theme::central_frame(&style));
        match self.mode {
            AppMode::Main => {
                self.handle_dropped_files(context);
                egui::SidePanel::right("restriction_panel")
                    .frame(theme::sidebar_frame(&style))
                    .resizable(true)
                    .default_width(360.0)
                    .min_width(280.0)
                    .show_animated(context, self.restriction_visible, |ui| {
                        self.restriction_panel_ui(ui);
                    });
                central.show(context, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            // Keep content clear of the floating scroll bar.
                            egui::Frame::new()
                                .inner_margin(Margin {
                                    right: 12,
                                    ..Margin::ZERO
                                })
                                .show(ui, |ui| self.main_analysis_ui(ui));
                        });
                });
                paint_drop_hint(context);
            }
            AppMode::Search => {
                central.show(context, |ui| self.search_window_ui(ui));
            }
            AppMode::Alignment => {
                central.show(context, |ui| self.alignment_window_ui(ui));
            }
        }
    }
}

/// Outline the window while a file is dragged over it.
fn paint_drop_hint(context: &egui::Context) {
    if context.input(|input| input.raw.hovered_files.is_empty()) {
        return;
    }
    let accent = Palette::for_dark_mode(context.style().visuals.dark_mode).accent;
    let painter = context.layer_painter(egui::LayerId::new(
        egui::Order::Foreground,
        egui::Id::new("file_drop_hint"),
    ));
    let rect = context.content_rect().shrink(4.0);
    painter.rect(
        rect,
        theme::RADIUS,
        accent.gamma_multiply(0.08),
        stroke(2.0, accent),
        egui::StrokeKind::Inside,
    );
    painter.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        "Drop a sequence file to load it",
        FontId::proportional(18.0),
        accent,
    );
}

/// Launch this executable again in a dedicated window mode, optionally
/// passing text through its standard input.
fn spawn_window(flag: &str, name: &str, input: Option<&str>) -> Result<Child, String> {
    let executable = std::env::current_exe()
        .map_err(|error| format!("Could not locate this application: {error}"))?;
    let mut command = Command::new(executable);
    command.arg(flag);
    if input.is_some() {
        command.stdin(Stdio::piped());
    }
    let mut child = command
        .spawn()
        .map_err(|error| format!("Could not open the {name} window: {error}"))?;
    if let Some(text) = input {
        let written = child
            .stdin
            .take()
            .ok_or_else(|| format!("The {name} window input pipe was unavailable."))
            .and_then(|mut pipe| {
                pipe.write_all(text.as_bytes()).map_err(|error| {
                    format!("Could not send sequence to the {name} window: {error}")
                })
            });
        if let Err(message) = written {
            let _ = child.kill();
            let _ = child.wait();
            return Err(message);
        }
    }
    Ok(child)
}

/// A read-only, selectable, scrollable monospace box.
fn output_box(ui: &mut egui::Ui, id_salt: &str, text: &str) {
    theme::surface_frame(ui).show(ui, |ui| {
        ui.set_min_width(ui.available_width());
        ScrollArea::vertical()
            .id_salt(id_salt)
            .max_height(OUTPUT_HEIGHT)
            .min_scrolled_height(OUTPUT_HEIGHT)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                // A `&str` buffer is immutable but still supports selection and copying.
                let mut text = text;
                ui.add(
                    TextEdit::multiline(&mut text)
                        .font(TextStyle::Monospace)
                        .frame(false)
                        .margin(Margin::ZERO)
                        .desired_rows(1)
                        .desired_width(f32::INFINITY),
                );
            });
    });
}

fn text_format(font: &FontId, color: Color32) -> TextFormat {
    TextFormat {
        font_id: font.clone(),
        color,
        ..Default::default()
    }
}

/// Append text, extending the previous section when the format is unchanged
/// so long sequences produce few sections.
fn append_run(job: &mut LayoutJob, text: &str, format: &TextFormat) {
    match job.sections.last_mut() {
        Some(last) if last.format == *format => {
            job.text.push_str(text);
            last.byte_range.end = job.text.len();
        }
        _ => job.append(text, 0.0, format.clone()),
    }
}

/// Append `bytes` as runs split wherever `flags` changes.
fn append_flagged(
    job: &mut LayoutJob,
    bytes: &[u8],
    flags: &[bool],
    plain: &TextFormat,
    flagged: &TextFormat,
) {
    let mut start = 0;
    while start < bytes.len() {
        let flag = flags[start];
        let end = flags[start..]
            .iter()
            .position(|&value| value != flag)
            .map_or(bytes.len(), |offset| start + offset);
        let text = std::str::from_utf8(&bytes[start..end]).unwrap_or("?");
        append_run(job, text, if flag { flagged } else { plain });
        start = end;
    }
}

fn origin_layout(
    sequence: &str,
    matches: &[usize],
    query_len: usize,
    palette: &Palette,
    font: FontId,
) -> LayoutJob {
    let bytes = sequence.as_bytes();
    let mut highlighted = vec![false; bytes.len()];
    for &start in matches {
        highlighted[start..(start + query_len).min(bytes.len())].fill(true);
    }

    let number = text_format(&font, palette.text_weak);
    let normal = text_format(&font, palette.text);
    let selected = TextFormat {
        background: palette.highlight_bg,
        ..text_format(&font, palette.highlight_fg)
    };
    let mut job = LayoutJob::default();
    for (line_index, line) in bytes.chunks(60).enumerate() {
        let line_start = line_index * 60;
        append_run(&mut job, &format!("{:>9} ", line_start + 1), &number);
        for (group_index, group) in line.chunks(10).enumerate() {
            let start = line_start + group_index * 10;
            if group_index > 0 {
                let joined = highlighted[start - 1] && highlighted[start];
                append_run(&mut job, " ", if joined { &selected } else { &normal });
            }
            append_flagged(
                &mut job,
                group,
                &highlighted[start..start + group.len()],
                &normal,
                &selected,
            );
        }
        append_run(&mut job, "\n", &normal);
    }
    job
}

/// Aligned rows with variable columns colored and a Clustal-style
/// conservation line (`*`) underneath.
fn alignment_layout(result: &AlignmentResult, palette: &Palette, font: FontId) -> LayoutJob {
    let label = text_format(&font, palette.text_weak);
    let normal = text_format(&font, palette.text);
    let mismatch = text_format(&font, palette.mismatch);
    let label_width = result.labels.iter().map(String::len).max().unwrap_or(4) + 2;
    let rows: Vec<&[u8]> = result.sequences.iter().map(String::as_bytes).collect();
    let width = rows.iter().map(|row| row.len()).max().unwrap_or(0);
    let variable: Vec<bool> = (0..width)
        .map(|column| {
            let first = rows.first().and_then(|row| row.get(column));
            rows.iter().any(|row| row.get(column) != first)
        })
        .collect();

    let mut job = LayoutJob::default();
    for (row_label, row) in result.labels.iter().zip(&rows) {
        append_run(&mut job, &format!("{row_label:<label_width$}"), &label);
        append_flagged(&mut job, row, &variable[..row.len()], &normal, &mismatch);
        append_run(&mut job, "\n", &normal);
    }
    let conservation: String = variable
        .iter()
        .map(|&is_variable| if is_variable { ' ' } else { '*' })
        .collect();
    append_run(&mut job, &" ".repeat(label_width), &label);
    append_run(&mut job, &conservation, &label);
    job
}

fn contains_ignore_ascii_case(haystack: &str, needle: &str) -> bool {
    needle.is_empty()
        || haystack
            .as_bytes()
            .windows(needle.len())
            .any(|window| window.eq_ignore_ascii_case(needle.as_bytes()))
}

fn group_thousands(value: usize) -> String {
    let digits = value.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_thousands() {
        assert_eq!(group_thousands(7), "7");
        assert_eq!(group_thousands(1_234), "1,234");
        assert_eq!(group_thousands(1_234_567), "1,234,567");
    }

    #[test]
    fn filters_case_insensitively() {
        assert!(contains_ignore_ascii_case("EcoRI", "ecor"));
        assert!(!contains_ignore_ascii_case("EcoRI", "bam"));
    }

    #[test]
    fn origin_layout_merges_runs() {
        let sequence = "A".repeat(120);
        let job = origin_layout(
            &sequence,
            &[5],
            3,
            Palette::for_dark_mode(false),
            FontId::monospace(13.0),
        );
        // Two lines of numbering plus bases, with one highlighted run.
        assert!(job.sections.len() < 10, "{} sections", job.sections.len());
        assert!(job.text.starts_with("        1 AAAAAAAAAA AAAAAAAAAA"));
    }
}
