#![cfg_attr(target_os = "windows", windows_subsystem = "windows")]

mod app;

use app::{AppMode, DnaAnalyzerApp};
use eframe::egui;
use std::io::Read;

fn main() -> eframe::Result {
    let arguments: Vec<_> = std::env::args_os().collect();
    let mode = if arguments
        .iter()
        .any(|argument| argument == "--search-window")
    {
        AppMode::Search
    } else if arguments
        .iter()
        .any(|argument| argument == "--alignment-window")
    {
        AppMode::Alignment
    } else {
        AppMode::Main
    };
    let initial_input = if mode == AppMode::Search {
        let mut sequence = String::new();
        if let Err(error) = std::io::stdin().read_to_string(&mut sequence) {
            eprintln!("Could not read sequence for search window: {error}");
        }
        Some(sequence)
    } else {
        None
    };
    let (title, inner_size, min_inner_size) = match mode {
        AppMode::Main => (
            "DNA/RNA Translator + Restriction Enzyme Analysis",
            [1500.0, 750.0],
            [900.0, 640.0],
        ),
        AppMode::Search => ("Sequence Search Highlight", [800.0, 600.0], [540.0, 360.0]),
        AppMode::Alignment => (
            "Multiple Sequence Alignment",
            [1000.0, 700.0],
            [680.0, 480.0],
        ),
    };
    let options = eframe::NativeOptions {
        // eframe 0.33.3 can crash during AppKit teardown on Intel Touch Bar
        // Macs when it persists native window state. This application does not
        // need that state, so disabling it avoids the duplicate KVO cleanup.
        persist_window: false,
        viewport: egui::ViewportBuilder::default()
            .with_title(title)
            .with_inner_size(inner_size)
            .with_min_inner_size(min_inner_size),
        ..Default::default()
    };

    eframe::run_native(
        title,
        options,
        Box::new(move |creation_context| {
            Ok(Box::new(DnaAnalyzerApp::new(
                creation_context,
                mode,
                initial_input,
            )))
        }),
    )
}
