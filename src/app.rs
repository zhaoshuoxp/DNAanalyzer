use dna_analyzer::alignment::{self, AlignmentResult};
use dna_analyzer::restriction::{self, RestrictionHit};
use dna_analyzer::{
    MoleculeKind, SequenceAnalysis, analyze_input, find_overlapping, format_origin_export,
    six_frame_translation,
};
use eframe::egui::{
    self, Color32, FontId, RichText, ScrollArea, TextEdit, TextFormat, TextStyle, text::LayoutJob,
};
use std::fs;
use std::io::Write;
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::time::{Duration, Instant};

struct RestrictionTask {
    receiver: Receiver<Vec<RestrictionHit>>,
    started: Instant,
    bases: usize,
}

struct AlignmentTask {
    receiver: Receiver<Result<AlignmentResult, String>>,
    started: Instant,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AppMode {
    Main,
    Search,
    Alignment,
}

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
    restriction_hits: Vec<RestrictionHit>,
    restriction_task: Option<RestrictionTask>,
    restriction_bases: Option<usize>,

    search_query: String,
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
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        creation_context
            .egui_ctx
            .style_mut(|style| style.spacing.item_spacing = egui::vec2(8.0, 8.0));

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
            restriction_hits: Vec::new(),
            restriction_task: None,
            restriction_bases: None,
            search_query: String::new(),
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
        self.search_matches = find_overlapping(&self.analysis.dna, &self.search_query);
        self.restriction_task = None;
        self.restriction_hits.clear();
        self.restriction_bases = None;
        if self.input.trim().is_empty() {
            self.status.clear();
            self.status_is_error = false;
        } else if self.analysis.dna.is_empty() {
            self.set_error("No valid DNA or RNA bases found.");
        } else {
            self.status = format!(
                "Parsed {} {} bases ({}).",
                self.analysis.dna.len(),
                self.analysis.kind.label(),
                self.analysis.input_format.label()
            );
            self.status_is_error = false;
        }
    }

    fn poll_tasks(&mut self, context: &egui::Context) {
        let restriction_update =
            self.restriction_task
                .as_ref()
                .and_then(|task| match task.receiver.try_recv() {
                    Ok(hits) => Some(Ok((hits, task.bases, task.started.elapsed()))),
                    Err(TryRecvError::Disconnected) => {
                        Some(Err("Restriction analysis worker stopped unexpectedly."))
                    }
                    Err(TryRecvError::Empty) => None,
                });
        if let Some(update) = restriction_update {
            self.restriction_task = None;
            match update {
                Ok((hits, bases, elapsed)) => {
                    self.status = format!(
                        "Found {} cutting enzymes in {:.1} ms.",
                        hits.len(),
                        elapsed.as_secs_f64() * 1_000.0
                    );
                    self.status_is_error = false;
                    self.restriction_hits = hits;
                    self.restriction_bases = Some(bases);
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
                    self.status = format!(
                        "Aligned {} sequences with {} in {:.2} s.",
                        alignment.sequences.len(),
                        alignment.engine,
                        elapsed.as_secs_f64()
                    );
                    self.status_is_error = false;
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
        let paths: Vec<_> = context.input(|input| {
            input
                .raw
                .dropped_files
                .iter()
                .filter_map(|file| file.path.clone())
                .collect()
        });
        if let Some(path) = paths.first() {
            self.load_file(path);
        }
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
                self.status = format!("Loaded {} bases.", self.analysis.dna.len());
                self.status_is_error = false;
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
                    self.status = format!("Saved numbered ORIGIN sequence to {}.", path.display());
                    self.status_is_error = false;
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
            Ok(()) => {
                self.status = "Opened NCBI BLAST in the default browser.".to_owned();
                self.status_is_error = false;
            }
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
            bases,
        });
        self.status = format!(
            "Scanning {bases} bases against {} enzymes…",
            restriction::enzymes().len()
        );
        self.status_is_error = false;
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
        self.status = "Running multiple sequence alignment…".to_owned();
        self.status_is_error = false;
    }

    fn open_search_window(&mut self) {
        if self.analysis.dna.is_empty() {
            self.set_error("Enter sequence first.");
            return;
        }
        let executable = match std::env::current_exe() {
            Ok(path) => path,
            Err(error) => {
                self.set_error(format!("Could not locate this application: {error}"));
                return;
            }
        };
        let mut child = match Command::new(executable)
            .arg("--search-window")
            .stdin(Stdio::piped())
            .spawn()
        {
            Ok(child) => child,
            Err(error) => {
                self.set_error(format!("Could not open the search window: {error}"));
                return;
            }
        };
        let write_result = child
            .stdin
            .take()
            .ok_or_else(|| "Search window input pipe was unavailable.".to_owned())
            .and_then(|mut input| {
                input
                    .write_all(self.analysis.display.as_bytes())
                    .map_err(|error| format!("Could not send sequence to search window: {error}"))
            });
        if let Err(message) = write_result {
            let _ = child.kill();
            let _ = child.wait();
            self.set_error(message);
            return;
        }
        self.child_processes.push(child);
        self.status = "Opened sequence search in a new window.".to_owned();
        self.status_is_error = false;
    }

    fn open_alignment_window(&mut self) {
        let executable = match std::env::current_exe() {
            Ok(path) => path,
            Err(error) => {
                self.set_error(format!("Could not locate this application: {error}"));
                return;
            }
        };
        match Command::new(executable).arg("--alignment-window").spawn() {
            Ok(child) => {
                self.child_processes.push(child);
                self.status = "Opened multiple alignment in a new window.".to_owned();
                self.status_is_error = false;
            }
            Err(error) => {
                self.set_error(format!("Could not open the alignment window: {error}"));
            }
        }
    }

    fn warning_line(&self, ui: &mut egui::Ui) {
        let mut warnings = Vec::new();
        if self.analysis.kind == MoleculeKind::Rna {
            warnings.push("RNA detected, using A:U base pairing.".to_owned());
        } else if self.analysis.kind == MoleculeKind::Mixed {
            warnings.push("Both T and U detected; U is normalized to T for analysis.".to_owned());
        }
        if !self.analysis.invalid_letters.is_empty() {
            let invalid: String = self.analysis.invalid_letters.iter().collect();
            warnings.push(format!("Invalid characters removed: {invalid}"));
        }

        let (message, error) = if warnings.is_empty() {
            (&self.status, self.status_is_error)
        } else {
            // This branch needs an owned joined value that lives through the label call.
            let joined = warnings.join(" | ");
            let color = ui.visuals().error_fg_color;
            ui.label(RichText::new(joined).color(color));
            return;
        };
        let color = if error {
            ui.visuals().error_fg_color
        } else {
            ui.visuals().weak_text_color()
        };
        ui.label(RichText::new(message).color(color));
    }

    fn main_analysis_ui(&mut self, ui: &mut egui::Ui) {
        self.warning_line(ui);
        ui.label("Enter DNA or RNA sequence (A/T/G/C/U):");
        let response = ScrollArea::vertical()
            .id_salt("main_sequence_input")
            .max_height(145.0)
            .min_scrolled_height(145.0)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                ui.add_sized(
                    [ui.available_width(), 145.0],
                    TextEdit::multiline(&mut self.input)
                        .font(TextStyle::Monospace)
                        .hint_text("Paste DNA, RNA, FASTA, or GenBank ORIGIN text here…")
                        .lock_focus(true),
                )
            })
            .inner;
        if response.changed() {
            self.refresh_analysis();
        }

        if full_width_button(ui, "Search Subsequence").clicked() {
            self.open_search_window();
        }

        let mut actions = [false; 6];
        ui.columns(6, |columns| {
            let labels = [
                "Clear",
                "Translate 6-Frames",
                "Show Restriction Enzymes",
                "NCBI BLAST",
                "Multiple Alignment",
                "Save as Fasta",
            ];
            for (index, column) in columns.iter_mut().enumerate() {
                actions[index] = column
                    .add_sized(
                        [column.available_width(), 28.0],
                        egui::Button::new(labels[index]),
                    )
                    .clicked();
            }
        });

        if actions[0] {
            self.input.clear();
            self.refresh_analysis();
        }
        if actions[1] {
            self.translation_visible = true;
        }
        if actions[2] {
            self.restriction_visible = true;
            self.start_restriction_analysis();
        }
        if actions[3] {
            self.open_blast();
        }
        if actions[4] {
            self.open_alignment_window();
        }
        if actions[5] {
            self.save_origin();
        }

        ui.add_space(5.0);
        ui.columns(3, |columns| {
            sequence_box(
                &mut columns[0],
                "Complement:",
                "Copy Complement",
                &self.analysis.complement,
            );
            sequence_box(
                &mut columns[1],
                "Reverse:",
                "Copy Reverse",
                &self.analysis.reverse,
            );
            sequence_box(
                &mut columns[2],
                "Reverse-Complement:",
                "Copy Reverse-Complement",
                &self.analysis.reverse_complement,
            );
        });

        if self.translation_visible {
            ui.add_space(5.0);
            for row in 0..2 {
                ui.columns(3, |columns| {
                    for (column, panel) in columns.iter_mut().enumerate() {
                        let index = row * 3 + column;
                        translation_box(panel, index + 1, &self.translations[index]);
                    }
                });
            }
        }
    }

    fn restriction_panel_ui(&mut self, ui: &mut egui::Ui) {
        ui.label(RichText::new("Restriction Enzyme Sites (1-based):").strong());
        ui.add(
            TextEdit::singleline(&mut self.restriction_filter)
                .hint_text("Filter enzymes (comma-separated)")
                .desired_width(f32::INFINITY),
        );
        if self.restriction_task.is_some() {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("Analyzing…");
            });
        }
        let needs_refresh = !self.analysis.dna.is_empty()
            && self.restriction_task.is_none()
            && self.restriction_bases != Some(self.analysis.dna.len());
        if needs_refresh {
            ui.label(
                RichText::new("Sequence changed; click Show Restriction Enzymes to refresh.")
                    .color(ui.visuals().warn_fg_color),
            );
        }

        let keywords: Vec<String> = self
            .restriction_filter
            .split(',')
            .map(|word| word.trim().to_ascii_uppercase())
            .filter(|word| !word.is_empty())
            .collect();
        let mut text = String::new();
        for hit in self.restriction_hits.iter().filter(|hit| {
            keywords.is_empty()
                || keywords.iter().any(|keyword| {
                    hit.name.to_ascii_uppercase().contains(keyword)
                        || hit.site.to_ascii_uppercase().contains(keyword)
                })
        }) {
            text.push_str(&format!(
                "{} ({}): {}\n",
                hit.name,
                hit.site,
                hit.positions
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        if text.is_empty() && self.restriction_task.is_none() {
            text = if self.analysis.dna.is_empty() {
                "Enter a sequence to analyze.".to_owned()
            } else if needs_refresh {
                "Restriction results are waiting to be refreshed.".to_owned()
            } else {
                "No cutting enzymes found.".to_owned()
            };
        }
        ui.add_sized(
            ui.available_size(),
            TextEdit::multiline(&mut text)
                .font(TextStyle::Monospace)
                .interactive(false),
        );
    }

    fn search_window_ui(&mut self, ui: &mut egui::Ui) {
        let response = ui.add(
            TextEdit::singleline(&mut self.search_query)
                .font(TextStyle::Monospace)
                .hint_text("Search sequence…")
                .desired_width(f32::INFINITY),
        );
        if response.changed() {
            self.search_matches = find_overlapping(&self.analysis.dna, &self.search_query);
        }
        let query_len = analyze_input(&self.search_query).dna.len();
        let layout = origin_layout(
            &self.analysis.display,
            &self.search_matches,
            query_len,
            ui.visuals().text_color(),
        );
        ScrollArea::both()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                ui.label(layout);
            });
    }

    fn alignment_window_ui(&mut self, ui: &mut egui::Ui) {
        self.warning_line(ui);
        ScrollArea::vertical()
            .id_salt("alignment_sequence_inputs")
            .max_height(310.0)
            .show(ui, |ui| {
                for sequence in &mut self.alignment_inputs {
                    ui.add_sized(
                        [ui.available_width(), 58.0],
                        TextEdit::multiline(sequence)
                            .font(TextStyle::Monospace)
                            .hint_text("Enter DNA/RNA or protein sequence…"),
                    );
                }
            });

        ui.columns(2, |columns| {
            if columns[0]
                .add_sized(
                    [columns[0].available_width(), 28.0],
                    egui::Button::new("Add Sequence"),
                )
                .clicked()
            {
                self.alignment_inputs.push(String::new());
            }
            if columns[1]
                .add_enabled_ui(self.alignment_inputs.len() > 2, |ui| {
                    ui.add_sized(
                        [ui.available_width(), 28.0],
                        egui::Button::new("Remove Sequence"),
                    )
                })
                .inner
                .clicked()
            {
                self.alignment_inputs.pop();
            }
        });

        if let Some(result) = &self.alignment_result {
            ui.label(
                RichText::new(&result.engine)
                    .small()
                    .color(ui.visuals().weak_text_color()),
            );
            ScrollArea::both()
                .id_salt("alignment_result")
                .max_height(260.0)
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(alignment_layout(result, ui.visuals().text_color()));
                });
        } else {
            egui::Frame::group(ui.style())
                .inner_margin(12.0)
                .show(ui, |ui| {
                    ui.set_min_height(220.0);
                    if self.alignment_task.is_some() {
                        ui.spinner();
                    }
                });
        }

        let running = self.alignment_task.is_some();
        if ui
            .add_enabled(
                !running,
                egui::Button::new(if running { "Aligning…" } else { "Align" })
                    .min_size(egui::vec2(ui.available_width(), 28.0)),
            )
            .clicked()
        {
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
        match self.mode {
            AppMode::Main => {
                self.handle_dropped_files(context);
                egui::SidePanel::right("restriction_panel")
                    .resizable(true)
                    .default_width(360.0)
                    .min_width(280.0)
                    .show_animated(context, self.restriction_visible, |ui| {
                        self.restriction_panel_ui(ui);
                    });
                egui::CentralPanel::default().show(context, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| self.main_analysis_ui(ui));
                });
            }
            AppMode::Search => {
                egui::CentralPanel::default().show(context, |ui| self.search_window_ui(ui));
            }
            AppMode::Alignment => {
                egui::CentralPanel::default().show(context, |ui| self.alignment_window_ui(ui));
            }
        }
    }
}

fn full_width_button(ui: &mut egui::Ui, text: &str) -> egui::Response {
    ui.add_sized([ui.available_width(), 28.0], egui::Button::new(text))
}

fn sequence_box(ui: &mut egui::Ui, title: &str, copy_label: &str, sequence: &str) {
    ui.label(title);
    let mut display = sequence.to_owned();
    ui.add_sized(
        [ui.available_width(), 105.0],
        TextEdit::multiline(&mut display)
            .font(TextStyle::Monospace)
            .interactive(false),
    );
    if full_width_button(ui, copy_label).clicked() {
        ui.ctx().copy_text(sequence.to_owned());
    }
}

fn translation_box(ui: &mut egui::Ui, frame: usize, protein: &str) {
    ui.label(format!("Frame {frame}:"));
    let mut display = protein.to_owned();
    ui.add_sized(
        [ui.available_width(), 105.0],
        TextEdit::multiline(&mut display)
            .font(TextStyle::Monospace)
            .interactive(false),
    );
}

fn origin_layout(
    sequence: &str,
    matches: &[usize],
    query_len: usize,
    text_color: Color32,
) -> LayoutJob {
    let mut highlighted = vec![false; sequence.len()];
    for &start in matches {
        for position in start..(start + query_len).min(highlighted.len()) {
            highlighted[position] = true;
        }
    }

    let normal = TextFormat {
        font_id: FontId::monospace(14.0),
        color: text_color,
        ..Default::default()
    };
    let selected = TextFormat {
        font_id: FontId::monospace(14.0),
        color: Color32::BLACK,
        background: Color32::YELLOW,
        ..Default::default()
    };
    let mut job = LayoutJob::default();
    job.wrap.max_width = f32::INFINITY;
    for line_start in (0..sequence.len()).step_by(60) {
        job.append(&format!("{:>9} ", line_start + 1), 0.0, normal.clone());
        let line_end = (line_start + 60).min(sequence.len());
        for position in line_start..line_end {
            if position > line_start && (position - line_start) % 10 == 0 {
                job.append(" ", 0.0, normal.clone());
            }
            job.append(
                &sequence[position..position + 1],
                0.0,
                if highlighted[position] {
                    selected.clone()
                } else {
                    normal.clone()
                },
            );
        }
        job.append("\n", 0.0, normal.clone());
    }
    job
}

fn alignment_layout(result: &AlignmentResult, text_color: Color32) -> LayoutJob {
    let normal = TextFormat {
        font_id: FontId::monospace(14.0),
        color: text_color,
        ..Default::default()
    };
    let mismatch = TextFormat {
        font_id: FontId::monospace(14.0),
        color: Color32::RED,
        ..Default::default()
    };
    let label_width = result.labels.iter().map(String::len).max().unwrap_or(4) + 2;
    let rows: Vec<&[u8]> = result.sequences.iter().map(String::as_bytes).collect();
    let mut job = LayoutJob::default();
    job.wrap.max_width = f32::INFINITY;
    for (row_index, row) in rows.iter().enumerate() {
        job.append(
            &format!("{:<label_width$}", result.labels[row_index]),
            0.0,
            normal.clone(),
        );
        for (column, base) in row.iter().enumerate() {
            let conserved = rows.iter().all(|candidate| candidate[column] == *base);
            job.append(
                std::str::from_utf8(&row[column..column + 1]).unwrap_or("?"),
                0.0,
                if conserved {
                    normal.clone()
                } else {
                    mismatch.clone()
                },
            );
        }
        job.append("\n", 0.0, normal.clone());
    }
    job
}
