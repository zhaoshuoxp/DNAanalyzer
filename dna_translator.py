#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QTextEdit, QPushButton, QLabel, QGridLayout, QLineEdit, QFileDialog, QDialog
)
from PyQt5.QtGui import QFont, QTextCursor, QTextCharFormat, QColor
from PyQt5.QtCore import Qt
from Bio.Restriction import AllEnzymes, RestrictionBatch
from Bio.Seq import Seq

# --- Genetic code ---
genetic_code = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

# --- Restriction enzymes batch ---
enzyme_batch = RestrictionBatch(AllEnzymes)

# --- Format sequence as ORIGIN / Fasta ---
def format_fasta_origin(seq, line_length=60, group=10):
    seq = seq.upper().replace("U","T")
    lines = []
    for i in range(0, len(seq), line_length):
        line_seq = seq[i:i+line_length]
        grouped = ' '.join([line_seq[j:j+group] for j in range(0,len(line_seq),group)])
        lines.append(f"{i+1:9} {grouped}")
    return "\n".join(lines)

# --- Sequence Search Highlight Dialog ---
class SequenceSearchDialog(QDialog):
    def __init__(self, seq, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sequence Search Highlight")
        self.resize(800, 600)
        layout = QVBoxLayout()
        
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search sequence...")
        self.search_input.textChanged.connect(self.highlight_sequence)
        layout.addWidget(self.search_input)
        
        self.text_area = QTextEdit()
        self.text_area.setFont(QFont("Courier", 12))
        layout.addWidget(self.text_area)
        self.setLayout(layout)
        
        # --- Clean sequence: remove invalid characters, U->T ---
        allowed = set("ACGTU")
        clean_seq = "".join([c.upper() for c in seq if c.upper() in allowed]).replace("U","T")
        self.seq = clean_seq
        
        # --- Format as ORIGIN / Fasta style ---
        self.formatted_seq = format_fasta_origin(self.seq)
        self.text_area.setText(self.formatted_seq)
        
    def highlight_sequence(self):
        search_text = self.search_input.text().upper()
        cursor = self.text_area.textCursor()
        cursor.select(QTextCursor.Document)
        cursor.setCharFormat(QTextCharFormat())
        cursor.clearSelection()
        
        if not search_text:
            return
        
        fmt = QTextCharFormat()
        fmt.setBackground(QColor("yellow"))
        fmt.setForeground(QColor("black"))
        
        text = self.text_area.toPlainText()
        start = 0
        while True:
            idx = text.find(search_text, start)
            if idx == -1:
                break
            cursor.setPosition(idx)
            cursor.movePosition(QTextCursor.Right, QTextCursor.KeepAnchor, len(search_text))
            cursor.mergeCharFormat(fmt)
            start = idx + len(search_text)
            
# --- Main GUI ---
class DNATranslator(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA/RNA Translator + Restriction Enzyme Analysis")
        self.setGeometry(100, 50, 1500, 750)
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        right_layout.setAlignment(Qt.AlignTop)
        
        # --- Warning label ---
        self.warning_label = QLabel("")
        self.warning_label.setStyleSheet("color: red;")
        left_layout.addWidget(self.warning_label)
        
        # --- Input sequence ---
        left_layout.addWidget(QLabel("Enter DNA or RNA sequence (A/T/G/C/U):"))
        self.input = QTextEdit()
        self.input.setFont(QFont("Courier", 12))
        left_layout.addWidget(self.input)
        
        # --- Sequence search highlighter ---
        self.seq_search_btn = QPushButton("Search Subsequence")
        self.seq_search_btn.clicked.connect(self.show_sequence_search_dialog)
        left_layout.addWidget(self.seq_search_btn)
        
        # --- Buttons ---
        btn_layout = QHBoxLayout()
        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(lambda: self.input.clear())
        self.translate_btn = QPushButton("Translate 6-Frames")
        self.translate_btn.clicked.connect(self.show_translation)
        self.re_btn = QPushButton("Show Restriction Enzymes")
        self.re_btn.clicked.connect(self.show_restriction_sites)
        self.save_btn = QPushButton("Save as Fasta")
        self.save_btn.clicked.connect(self.save_fasta)
        btn_layout.addWidget(self.clear_btn)
        btn_layout.addWidget(self.translate_btn)
        btn_layout.addWidget(self.re_btn)
        btn_layout.addWidget(self.save_btn)
        left_layout.addLayout(btn_layout)
        
        # --- Complement / Reverse / Reverse-Complement ---
        seq_layout = QGridLayout()
        self.comp_label = QTextEdit(); self.comp_label.setReadOnly(True); self.comp_label.setFont(QFont("Courier",12))
        self.rev_label = QTextEdit(); self.rev_label.setReadOnly(True); self.rev_label.setFont(QFont("Courier",12))
        self.revc_label = QTextEdit(); self.revc_label.setReadOnly(True); self.revc_label.setFont(QFont("Courier",12))
        
        self.comp_copy = QPushButton("Copy Complement"); self.comp_copy.clicked.connect(lambda: self.copy_text(self.comp_label.toPlainText()))
        self.rev_copy = QPushButton("Copy Reverse"); self.rev_copy.clicked.connect(lambda: self.copy_text(self.rev_label.toPlainText()))
        self.revc_copy = QPushButton("Copy Reverse-Complement"); self.revc_copy.clicked.connect(lambda: self.copy_text(self.revc_label.toPlainText()))
        
        seq_layout.addWidget(QLabel("Complement:"),0,0); seq_layout.addWidget(self.comp_label,1,0); seq_layout.addWidget(self.comp_copy,2,0)
        seq_layout.addWidget(QLabel("Reverse:"),0,1); seq_layout.addWidget(self.rev_label,1,1); seq_layout.addWidget(self.rev_copy,2,1)
        seq_layout.addWidget(QLabel("Reverse-Complement:"),0,2); seq_layout.addWidget(self.revc_label,1,2); seq_layout.addWidget(self.revc_copy,2,2)
        left_layout.addLayout(seq_layout)
        
        # --- 6-frame translation ---
        self.translation_container = QWidget()
        self.translation_container.hide()
        trans_layout = QGridLayout()
        self.translation_labels = []
        for i in range(6):
            txt = QTextEdit(); txt.setReadOnly(True); txt.setFont(QFont("Courier",12))
            self.translation_labels.append(txt)
            trans_layout.addWidget(QLabel(f"Frame {i+1}:"), i//3, (i%3)*2)
            trans_layout.addWidget(txt, i//3, (i%3)*2+1)
        self.translation_container.setLayout(trans_layout)
        left_layout.addWidget(self.translation_container)
        
        # --- Restriction enzyme analysis ---
        self.re_label_title = QLabel("Restriction Enzyme Sites (1-based):")
        self.re_label_title.hide()
        self.re_search = QLineEdit()
        self.re_search.setPlaceholderText("Filter enzymes (e.g., EcoRI)")
        self.re_search.hide()
        self.re_search.textChanged.connect(self.update_restriction_sites)
        self.re_enzyme_label = QTextEdit()
        self.re_enzyme_label.setReadOnly(True)
        self.re_enzyme_label.setFont(QFont("Courier",12))
        self.re_enzyme_label.hide()
        right_layout.addWidget(self.re_label_title)
        right_layout.addWidget(self.re_search)
        right_layout.addWidget(self.re_enzyme_label)
        
        # --- Combine layouts ---
        main_layout.addLayout(left_layout,3)
        main_layout.addLayout(right_layout,1)
        self.setLayout(main_layout)
        
        # --- Input change ---
        self.input.textChanged.connect(self.update_sequences)
        
    # -------------------
    def update_sequences(self):
        seq = self.input.toPlainText().upper()
        has_U = 'U' in seq
        allowed = set("ACGTU")
        invalid_chars = set(seq) - allowed
        messages = []
        if has_U:
            messages.append("RNA detected, using A:U base pairing.")
        if invalid_chars:
            messages.append(f"Invalid characters removed: {''.join(invalid_chars)}")
        self.warning_label.setText(" | ".join(messages))
        clean_seq = "".join([c for c in seq if c in allowed]).replace('U','T')
        self.comp_label.setText(self.complement(clean_seq, has_U))
        self.rev_label.setText(self.reverse(clean_seq))
        self.revc_label.setText(self.reverse_complement(clean_seq, has_U))
        
    def complement(self, seq, is_rna=False):
        table = str.maketrans("ACGT", "UGCA" if is_rna else "TGCA")
        return seq.translate(table)
    
    def reverse(self, seq): return seq[::-1]
    def reverse_complement(self, seq, is_rna=False): return self.complement(seq,is_rna)[::-1]
        
    # -------------------
    # 6-frame translation
    def show_translation(self):
        self.translation_container.show()
        seq = self.input.toPlainText().upper()
        clean_seq = "".join([c for c in seq if c in "ACGTU"]).replace("U","T")
        revc = self.reverse_complement(clean_seq)
        frames = [self.translate(clean_seq,f) for f in range(3)] + [self.translate(revc,f) for f in range(3)]
        for i,lbl in enumerate(self.translation_labels): lbl.setText(frames[i])
        
    def translate(self, seq, frame):
        seq = seq[frame:]
        aa = []
        for i in range(0,len(seq)-2,3):
            codon = seq[i:i+3]; aa.append(genetic_code.get(codon,"X"))
        return "".join(aa)
    
    # -------------------
    # Copy to clipboard
    def copy_text(self,text): QApplication.clipboard().setText(text)
        
    # -------------------
    # Restriction enzyme analysis
    def show_restriction_sites(self):
        self.re_label_title.show()
        self.re_search.show()
        self.re_enzyme_label.show()
        self.update_restriction_sites()
        
    def update_restriction_sites(self):
        seq = "".join([c for c in self.input.toPlainText().upper() if c in "ACGTU"]).replace("U","T")
        seq_obj = Seq(seq)
        keyword = self.re_search.text().upper().strip()
        results = []
        for enzyme in enzyme_batch:
            if keyword and keyword not in enzyme.__name__.upper():
                continue
            positions = [pos+1 for pos in enzyme.search(seq_obj)]
            if positions:
                results.append(f"{enzyme.__name__} ({enzyme.site}): {', '.join(map(str, positions))}")
        self.re_enzyme_label.setText("\n".join(results))
        
    # -------------------
    # Save fasta
    def save_fasta(self):
        seq = "".join([c for c in self.input.toPlainText().upper() if c in "ACGTU"]).replace("U","T")
        formatted = format_fasta_origin(seq)
        filename, _ = QFileDialog.getSaveFileName(self, "Save Fasta", "", "Fasta Files (*.fasta *.fa);;All Files (*)")
        if filename:
            with open(filename, 'w') as f:
                f.write(">sequence_1\n")
                f.write(formatted + "\n//")
                
    # -------------------
    # Show sequence search dialog
    def show_sequence_search_dialog(self):
        seq = self.input.toPlainText()
        if not seq.strip():
            return
        dialog = SequenceSearchDialog(seq, self)
        dialog.exec_()
        
if __name__=="__main__":
    app = QApplication(sys.argv)
    win = DNATranslator()
    win.show()
    sys.exit(app.exec_())
    