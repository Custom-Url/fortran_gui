import sys
import json
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QMenu, QTextEdit, QMessageBox, QTabWidget, QWidget, QVBoxLayout, QLabel, QPushButton
)
from PyQt6.QtGui import QAction

from material_manager import MaterialManager

class FileManagerApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("FORTRAN Project GUI")
        self.resize(800, 600)

        self.recent_file_store = "recent_files.json"
        self.recent_files = self.load_recent_files()

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Add tabs
        self.tabs.addTab(self.material_tab(), "Material Properties")
        self.tabs.addTab(self.loading_tab(), "Loading Conditions")
        self.tabs.addTab(self.make_output_tab(), "Results")

        # Build menu
        self.menu = self.menuBar()

        file_menu = self.menu.addMenu("File")

        open_action = QAction("Open File", self)
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        save_proj_action = QAction("Save Project", self)
        save_proj_action.triggered.connect(self.save_project)
        file_menu.addAction(save_proj_action)

        load_proj_action = QAction("Load Project", self)
        load_proj_action.triggered.connect(self.load_project)
        file_menu.addAction(load_proj_action)

        # Recent files submenu
        self.recent_menu = QMenu("Open Recent", self)
        file_menu.addMenu(self.recent_menu)
        self.update_recent_menu()

    # ------------------------
    # Tab Layouts
    # ------------------------
    def material_tab(self):
        """Tab for file operations"""
        tab = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("File operations here"))
        layout.addWidget(QPushButton("New Material File"))

        btn_open_mat = QPushButton("Open Material File")
        btn_open_mat.clicked.connect(self.open_material_file)
        layout.addWidget(btn_open_mat)

        btn_save_mat = QPushButton("Save Material File")
        btn_save_mat.clicked.connect(self.save_material_file)
        layout.addWidget(btn_save_mat)

        self.text_area = QTextEdit()
        self.text_area.setReadOnly(True)
        self.text_area.setPlaceholderText("Material file contents will appear here...")
        layout.addWidget(self.text_area)

        tab.setLayout(layout)
        return tab

    def loading_tab(self):
        """Tab for project settings"""
        tab = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Project settings here"))
        layout.addWidget(QTextEdit("Parameters..."))

        tab.setLayout(layout)
        return tab

    def make_output_tab(self):
        """Tab for output/logging"""
        tab = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Output log"))
        layout.addWidget(QTextEdit())

        tab.setLayout(layout)
        return tab


    # ------------------------
    # File handling
    # ------------------------
    def open_file(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open File")
        if filepath:
            #self.text_area.append(f"Opened: {filepath}")
            self.add_to_recent(filepath)

    def save_project(self):
        project = {"files": self.recent_files}
        filepath, _ = QFileDialog.getSaveFileName(self, "Save Project", filter="Project Files (*.proj)")
        if filepath:
            with open(filepath, "w") as f:
                json.dump(project, f)
            QMessageBox.information(self, "Saved", f"Project saved to {filepath}")

    def load_project(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Load Project", filter="Project Files (*.proj)")
        if filepath:
            with open(filepath, "r") as f:
                project = json.load(f)
            self.recent_files = project.get("files", [])
            self.update_recent_menu()
            QMessageBox.information(self, "Loaded", f"Loaded project {filepath}")

    def open_material_file(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open Material XML", filter="XML Files (*.xml)")
        if filepath:
            self.material_manager = MaterialManager(filepath)
            self.material_manager.pretty_print()  # just for testing
            self.text_area.append(f"Loaded material file: {filepath}")

    def save_material_file(self):
        if hasattr(self, "material_manager"):
            self.material_manager.save_file()
            self.text_area.append(f"Saved material file: {self.material_manager.filepath}")
    # ------------------------
    # Recent files logic
    # ------------------------
    def add_to_recent(self, filepath):
        if filepath not in self.recent_files:
            self.recent_files.insert(0, filepath)
        self.recent_files = self.recent_files[:5]
        self.save_recent_files()
        self.update_recent_menu()

    def update_recent_menu(self):
        self.recent_menu.clear()
        for f in self.recent_files:
            action = QAction(f, self)
            action.triggered.connect(lambda checked, path=f: self.open_recent(path))
            self.recent_menu.addAction(action)

    def open_recent(self, filepath):
        if os.path.exists(filepath):
            self.text_area.append(f"Opened recent: {filepath}")
        else:
            QMessageBox.warning(self, "File not found", f"{filepath} does not exist!")

    def save_recent_files(self):
        with open(self.recent_file_store, "w") as f:
            json.dump(self.recent_files, f)

    def load_recent_files(self):
        if os.path.exists(self.recent_file_store):
            with open(self.recent_file_store, "r") as f:
                return json.load(f)
        return []


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FileManagerApp()
    window.show()
    sys.exit(app.exec())


