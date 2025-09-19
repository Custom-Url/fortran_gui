import sys
import json
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QMenu, QTextEdit, QMessageBox, QTabWidget, QWidget, QVBoxLayout, QLabel, QPushButton, QListWidget, QComboBox, QGroupBox, QFormLayout, QLineEdit
)
from PyQt6.QtGui import QAction

import material_manager as mm

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
        self.tabs.addTab(self.simulation_tab(), "Simulation")

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
        """Tab for selecting up to three materials"""
        tab = QWidget()
        layout = QVBoxLayout()

        def create_material_section(title):
            group = QGroupBox(title)
            group_layout = QVBoxLayout()

            combo = QComboBox()
            group_layout.addWidget(combo)

            form_layout = QFormLayout()
            name_field = QLineEdit()
            name_field.setReadOnly(True)
            lambda_field = QLineEdit()
            lambda_field.setReadOnly(True)
            mu_field = QLineEdit()
            mu_field.setReadOnly(True)
            form_layout.addRow("Name:", name_field)
            form_layout.addRow("Lambda:", lambda_field)
            form_layout.addRow("Mu:", mu_field)
            group_layout.addLayout(form_layout)

            group.setLayout(group_layout)
            return group, combo, {"name": name_field, "lambda": lambda_field, "mu": mu_field}

        # Build 3 sections
        self.material_sections = []
        for i in range(1, 4):
            section, combo, fields = create_material_section(f"Material {i}")
            layout.addWidget(section)
            self.material_sections.append((combo, fields))

        # Button to open Material Manager window
        btn_manage = QPushButton("Manage Materials")
        btn_manage.clicked.connect(self.open_material_manager)
        layout.addWidget(btn_manage)

        # Load materials into dropdowns
        self.refresh_materials()

        tab.setLayout(layout)
        return tab


    def open_material_manager(self):
        dialog = mm.MaterialManager(self)   # open the dialog
        dialog.exec()                       # wait until user closes it
        self.refresh_materials()            # reload list after closing


    def loading_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Project settings here"))
        layout.addWidget(QTextEdit("Parameters..."))

        tab.setLayout(layout)
        return tab

    def simulation_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Output log"))
        layout.addWidget(QTextEdit())

        btn_open_mat = QPushButton("Execute")
        btn_open_mat.clicked.connect(self.open_material_file)
        layout.addWidget(btn_open_mat)

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

    # -------------------------
    # GUI Helper Functions
    # -------------------------
    def refresh_materials(self):
        """Reload all dropdowns from material_manager"""
        materials = mm.load_materials()

        for combo, fields in self.material_sections:
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("Select material...")
            for mat in materials:
                combo.addItem(mat["Name"])
            combo.blockSignals(False)

            combo.currentIndexChanged.connect(
                lambda idx, c=combo, f=fields: self.fill_fields_from_selection(c, f)
            )

    def fill_fields_from_selection(self, combo, fields):
        name = combo.currentText()
        if name == "Select material...":
            for field in fields.values():
                field.clear()
            return

        mat = mm.get_material(name)
        if mat:
            fields["name"].setText(mat["Name"])
            fields["lambda"].setText(mat["Lambda"])
            fields["mu"].setText(mat["Mu"])


    def gui_create_material(self, fields):
        new_data = {
            "Name": fields["name"].text(),
            "Lambda": fields["lambda"].text(),
            "Mu": fields["mu"].text(),
        }
        mm.create_material(new_data)
        QMessageBox.information(None, "Material Added", f"Material '{new_data['Name']}' added.")
        self.refresh_materials()


    def gui_edit_material(self, fields):
        updated_data = {
            "Name": fields["name"].text(),
            "Lambda": fields["lambda"].text(),
            "Mu": fields["mu"].text(),
        }
        mm.edit_material(updated_data)
        QMessageBox.information(None, "Material Updated", f"Material '{updated_data['Name']}' updated.")
        self.refresh_materials()


    def gui_delete_material(self, combo):
        name_to_delete = combo.currentText()
        if name_to_delete == "Select material...":
            return
        mm.delete_material(name_to_delete)
        QMessageBox.information(None, "Material Deleted", f"Material '{name_to_delete}' deleted.")
        self.refresh_materials()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FileManagerApp()
    window.show()
    sys.exit(app.exec())


