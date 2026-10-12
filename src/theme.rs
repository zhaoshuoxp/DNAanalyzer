//! Light and dark visual themes. egui switches between them automatically
//! when the operating system appearance changes.

use eframe::egui::{
    self, Color32, CornerRadius, FontFamily, FontId, Frame, Margin, Shadow, Stroke, TextStyle,
    Theme, ThemePreference, Visuals,
};

pub const RADIUS: u8 = 6;

/// Semantic colors shared by the egui style and custom-painted text.
#[derive(Clone, Copy)]
pub struct Palette {
    pub background: Color32,
    pub sidebar: Color32,
    pub surface: Color32,
    pub control: Color32,
    pub control_hover: Color32,
    pub control_active: Color32,
    pub border: Color32,
    pub border_strong: Color32,
    pub text: Color32,
    pub text_strong: Color32,
    pub text_weak: Color32,
    pub accent: Color32,
    pub accent_hover: Color32,
    pub on_accent: Color32,
    pub warn: Color32,
    pub error: Color32,
    pub highlight_bg: Color32,
    pub highlight_fg: Color32,
    pub mismatch: Color32,
}

const LIGHT: Palette = Palette {
    background: Color32::from_rgb(0xF4, 0xF6, 0xF8),
    sidebar: Color32::from_rgb(0xEC, 0xEF, 0xF3),
    surface: Color32::from_rgb(0xFF, 0xFF, 0xFF),
    control: Color32::from_rgb(0xFF, 0xFF, 0xFF),
    control_hover: Color32::from_rgb(0xEE, 0xF1, 0xF4),
    control_active: Color32::from_rgb(0xE2, 0xE6, 0xEB),
    border: Color32::from_rgb(0xD6, 0xDA, 0xE0),
    border_strong: Color32::from_rgb(0xB4, 0xBB, 0xC5),
    text: Color32::from_rgb(0x2A, 0x2F, 0x36),
    text_strong: Color32::from_rgb(0x11, 0x14, 0x18),
    text_weak: Color32::from_rgb(0x6A, 0x73, 0x7F),
    accent: Color32::from_rgb(0x0B, 0x7A, 0x87),
    accent_hover: Color32::from_rgb(0x09, 0x67, 0x72),
    on_accent: Color32::WHITE,
    warn: Color32::from_rgb(0xA8, 0x5D, 0x00),
    error: Color32::from_rgb(0xC4, 0x2B, 0x2B),
    highlight_bg: Color32::from_rgb(0xFF, 0xDF, 0x6E),
    highlight_fg: Color32::from_rgb(0x1A, 0x16, 0x00),
    mismatch: Color32::from_rgb(0xCF, 0x2E, 0x3A),
};

const DARK: Palette = Palette {
    background: Color32::from_rgb(0x18, 0x1B, 0x20),
    sidebar: Color32::from_rgb(0x14, 0x16, 0x1A),
    surface: Color32::from_rgb(0x0F, 0x11, 0x14),
    control: Color32::from_rgb(0x25, 0x29, 0x30),
    control_hover: Color32::from_rgb(0x2E, 0x33, 0x3B),
    control_active: Color32::from_rgb(0x38, 0x3E, 0x47),
    border: Color32::from_rgb(0x2F, 0x34, 0x3C),
    border_strong: Color32::from_rgb(0x4A, 0x51, 0x5B),
    text: Color32::from_rgb(0xD9, 0xDD, 0xE2),
    text_strong: Color32::from_rgb(0xF5, 0xF7, 0xF9),
    text_weak: Color32::from_rgb(0x8E, 0x96, 0xA0),
    accent: Color32::from_rgb(0x3C, 0xC0, 0xCC),
    accent_hover: Color32::from_rgb(0x66, 0xD2, 0xDB),
    on_accent: Color32::from_rgb(0x06, 0x1B, 0x1E),
    warn: Color32::from_rgb(0xE9, 0xB0, 0x4D),
    error: Color32::from_rgb(0xFF, 0x7B, 0x72),
    highlight_bg: Color32::from_rgb(0xE3, 0xB3, 0x41),
    highlight_fg: Color32::from_rgb(0x14, 0x10, 0x00),
    mismatch: Color32::from_rgb(0xFF, 0x80, 0x78),
};

impl Palette {
    pub fn for_dark_mode(dark_mode: bool) -> &'static Self {
        if dark_mode { &DARK } else { &LIGHT }
    }

    pub fn of(ui: &egui::Ui) -> &'static Self {
        Self::for_dark_mode(ui.visuals().dark_mode)
    }
}

/// Install both themes and follow the system appearance.
///
/// `DNA_ANALYZER_THEME=light|dark` forces one appearance.
pub fn install(context: &egui::Context) {
    install_fonts(context);
    context.set_style_of(Theme::Light, style(Visuals::light(), &LIGHT));
    context.set_style_of(Theme::Dark, style(Visuals::dark(), &DARK));
    let preference = match std::env::var("DNA_ANALYZER_THEME").as_deref() {
        Ok("light") => ThemePreference::Light,
        Ok("dark") => ThemePreference::Dark,
        _ => ThemePreference::System,
    };
    context.set_theme(preference);
}

/// Font files tried in order; the first readable, parseable one wins.
type FontCandidates = &'static [(&'static str, u32)];

#[cfg(target_os = "macos")]
const SYSTEM_FONTS: [FontCandidates; 3] = [
    &[
        ("/System/Library/Fonts/SFNS.ttf", 0),
        ("/System/Library/Fonts/HelveticaNeue.ttc", 0),
    ],
    &[("/System/Library/Fonts/Menlo.ttc", 0)],
    &[
        ("/System/Library/Fonts/Hiragino Sans GB.ttc", 0),
        ("/System/Library/Fonts/STHeiti Light.ttc", 0),
    ],
];

#[cfg(target_os = "windows")]
const SYSTEM_FONTS: [FontCandidates; 3] = [
    &[
        ("C:\\Windows\\Fonts\\segoeui.ttf", 0),
        ("C:\\Windows\\Fonts\\arial.ttf", 0),
    ],
    &[
        ("C:\\Windows\\Fonts\\consola.ttf", 0),
        ("C:\\Windows\\Fonts\\cour.ttf", 0),
    ],
    &[
        ("C:\\Windows\\Fonts\\msyh.ttc", 0),
        ("C:\\Windows\\Fonts\\simsun.ttc", 0),
    ],
];

#[cfg(not(any(target_os = "macos", target_os = "windows")))]
const SYSTEM_FONTS: [FontCandidates; 3] = [
    &[
        ("/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf", 0),
        ("/usr/share/fonts/noto/NotoSans-Regular.ttf", 0),
        ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 0),
        ("/usr/share/fonts/TTF/DejaVuSans.ttf", 0),
    ],
    &[
        ("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 0),
        ("/usr/share/fonts/TTF/DejaVuSansMono.ttf", 0),
        ("/usr/share/fonts/truetype/noto/NotoSansMono-Regular.ttf", 0),
    ],
    &[
        ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", 0),
        ("/usr/share/fonts/noto-cjk/NotoSansCJK-Regular.ttc", 0),
        ("/usr/share/fonts/truetype/wqy/wqy-microhei.ttc", 0),
    ],
];

/// Prefer the platform's regular-weight UI and monospace fonts over egui's
/// light default, with a CJK fallback for file paths and names. egui's
/// bundled fonts stay at the end of each family for anything else.
fn install_fonts(context: &egui::Context) {
    let [proportional, monospace, cjk] = SYSTEM_FONTS.map(load_font);
    let mut fonts = egui::FontDefinitions::default();
    for (family, primary, name) in [
        (FontFamily::Proportional, proportional, "system-ui"),
        (FontFamily::Monospace, monospace, "system-mono"),
    ] {
        let mut names = Vec::new();
        if let Some(data) = primary {
            fonts.font_data.insert(name.to_owned(), data.into());
            names.push(name.to_owned());
        }
        if cjk.is_some() {
            names.push("system-cjk".to_owned());
        }
        if let Some(list) = fonts.families.get_mut(&family) {
            list.splice(0..0, names);
        }
    }
    if let Some(data) = cjk {
        // CJK faces sit high next to Latin UI fonts; lower them slightly.
        let data = data.tweak(egui::FontTweak {
            y_offset_factor: 0.15,
            ..Default::default()
        });
        fonts.font_data.insert("system-cjk".to_owned(), data.into());
    }
    context.set_fonts(fonts);
}

fn load_font(candidates: FontCandidates) -> Option<egui::FontData> {
    candidates.iter().find_map(|&(path, index)| {
        let bytes = std::fs::read(path).ok()?;
        // egui panics on fonts it cannot parse, so validate them here first.
        let font = ab_glyph::FontRef::try_from_slice_and_index(&bytes, index).ok()?;
        ab_glyph::Font::units_per_em(&font)?;
        let mut data = egui::FontData::from_owned(bytes);
        data.index = index;
        Some(data)
    })
}

fn style(mut visuals: Visuals, palette: &Palette) -> egui::Style {
    let radius = CornerRadius::same(RADIUS);
    let border = stroke(1.0, palette.border);

    visuals.panel_fill = palette.background;
    visuals.window_fill = palette.background;
    visuals.window_stroke = border;
    visuals.window_corner_radius = CornerRadius::same(10);
    visuals.menu_corner_radius = radius;
    visuals.extreme_bg_color = palette.surface;
    visuals.text_edit_bg_color = Some(palette.surface);
    visuals.faint_bg_color = palette.control_hover;
    visuals.code_bg_color = palette.control_hover;
    visuals.weak_text_color = Some(palette.text_weak);
    visuals.hyperlink_color = palette.accent;
    visuals.warn_fg_color = palette.warn;
    visuals.error_fg_color = palette.error;
    visuals.selection.bg_fill = palette.accent.gamma_multiply(0.28);
    visuals.selection.stroke = stroke(1.5, palette.accent);
    visuals.text_cursor.stroke = stroke(2.0, palette.accent);
    visuals.popup_shadow = Shadow {
        offset: [0, 4],
        blur: 12,
        spread: 0,
        color: Color32::from_black_alpha(if visuals.dark_mode { 90 } else { 28 }),
    };
    visuals.window_shadow = visuals.popup_shadow;

    let widgets = &mut visuals.widgets;
    widgets.noninteractive.bg_fill = palette.background;
    widgets.noninteractive.weak_bg_fill = palette.background;
    widgets.noninteractive.bg_stroke = border;
    widgets.noninteractive.fg_stroke = stroke(1.0, palette.text);
    for (state, fill, border_color) in [
        (&mut widgets.inactive, palette.control, palette.border),
        (
            &mut widgets.hovered,
            palette.control_hover,
            palette.border_strong,
        ),
        (&mut widgets.active, palette.control_active, palette.accent),
        (
            &mut widgets.open,
            palette.control_active,
            palette.border_strong,
        ),
    ] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.bg_stroke = stroke(1.0, border_color);
        state.fg_stroke = stroke(1.0, palette.text_strong);
        state.expansion = 0.0;
    }
    for state in [
        &mut widgets.noninteractive,
        &mut widgets.inactive,
        &mut widgets.hovered,
        &mut widgets.active,
        &mut widgets.open,
    ] {
        state.corner_radius = radius;
    }

    let mut style = egui::Style {
        visuals,
        ..Default::default()
    };
    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(10.0, 5.0);
    style.spacing.interact_size.y = 26.0;
    style.text_styles = [
        (
            TextStyle::Small,
            FontId::new(11.0, FontFamily::Proportional),
        ),
        (TextStyle::Body, FontId::new(13.5, FontFamily::Proportional)),
        (
            TextStyle::Button,
            FontId::new(13.5, FontFamily::Proportional),
        ),
        (
            TextStyle::Heading,
            FontId::new(17.0, FontFamily::Proportional),
        ),
        (
            TextStyle::Monospace,
            FontId::new(13.0, FontFamily::Monospace),
        ),
    ]
    .into();
    style
}

/// A rounded, bordered surface used behind read-only sequence output.
pub fn surface_frame(ui: &egui::Ui) -> Frame {
    let palette = Palette::of(ui);
    Frame::new()
        .fill(palette.surface)
        .stroke(stroke(1.0, palette.border))
        .corner_radius(RADIUS)
        .inner_margin(Margin::symmetric(8, 6))
}

pub fn central_frame(style: &egui::Style) -> Frame {
    Frame::central_panel(style).inner_margin(Margin::same(14))
}

pub fn sidebar_frame(style: &egui::Style) -> Frame {
    Frame::side_top_panel(style)
        .fill(Palette::for_dark_mode(style.visuals.dark_mode).sidebar)
        .inner_margin(Margin::same(14))
}

/// A full-width button. `primary` buttons are accent-filled.
pub fn wide_button(ui: &mut egui::Ui, enabled: bool, primary: bool, text: &str) -> egui::Response {
    let size = [ui.available_width(), 28.0];
    if !primary {
        return ui
            .add_enabled_ui(enabled, |ui| ui.add_sized(size, egui::Button::new(text)))
            .inner;
    }
    let palette = Palette::of(ui);
    ui.scope(|ui| {
        let widgets = &mut ui.visuals_mut().widgets;
        for (state, fill) in [
            (&mut widgets.inactive, palette.accent),
            (&mut widgets.hovered, palette.accent_hover),
            (&mut widgets.active, palette.accent_hover),
        ] {
            state.weak_bg_fill = fill;
            state.bg_stroke = Stroke::NONE;
            state.fg_stroke = stroke(1.0, palette.on_accent);
        }
        ui.add_enabled_ui(enabled, |ui| ui.add_sized(size, egui::Button::new(text)))
            .inner
    })
    .inner
}

/// A muted caption above a field.
pub fn caption(ui: &mut egui::Ui, text: &str) {
    let color = Palette::of(ui).text_weak;
    ui.label(egui::RichText::new(text).size(12.5).color(color));
}

pub fn stroke(width: f32, color: Color32) -> Stroke {
    Stroke::new(width, color)
}
