import json
from PySide6 import QtCore, QtGui, QtWidgets
import plotly.graph_objects as go
import numpy as np
from eval import *

def isprimitive(data):
    # print("testing ", data) 
    if data is None or isinstance(data, str) or isinstance(data,int) or isinstance(data, float):
        return True
    if isinstance(data, bool):
        return True
    return False

class JsonViewerWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("JSON Viewer")

        # Create layout
        layout = QtWidgets.QVBoxLayout(self)

        # Create a QTreeView to display JSON
        self.tree_view = QtWidgets.QTreeView(self)
        self.tree_view.setAlternatingRowColors(True)
        self.tree_view.setRootIsDecorated(True)  # Display arrows for expanding/collapsing
        self.tree_view.setUniformRowHeights(True)
        layout.addWidget(self.tree_view)

        # Create a QPushButton to open a JSON file
        # self.open_button = QtWidgets.QPushButton("Open JSON File", self)
        # self.open_button.clicked.connect(self.open_json_file)
        # layout.addWidget(self.open_button)

        self.plot_container = QtWidgets.QWidget(self)
        layout.addWidget(self.plot_container)
        plot_layout = QtWidgets.QHBoxLayout(self.plot_container)
        self.plot_container.setLayout(plot_layout)
        self.plot_monom_button = QtWidgets.QPushButton("Plot monomial", self)
        self.plot_monom_button.clicked.connect(self.plot_selected_as_monomial)
        plot_layout.addWidget(self.plot_monom_button)
        
        self.plot_cheb_button = QtWidgets.QPushButton("Plot chebyshev", self)
        self.plot_cheb_button.clicked.connect(self.plot_selected_as_chebyshev)
        plot_layout.addWidget(self.plot_cheb_button)
        
        # Create a label for showing error messages
        self.error_label = QtWidgets.QLabel("", self)
        self.error_label.setStyleSheet("color: red;")
        layout.addWidget(self.error_label)

        self.holdToggle = QtWidgets.QCheckBox("Hold figure", self)
        layout.addWidget(self.holdToggle)
        def holdHandler(state):
            self.holdFigure = QtCore.Qt.CheckState(state) == QtCore.Qt.Checked
            print(self.holdFigure)
        self.holdToggle.stateChanged.connect(holdHandler)

        # Set the layout
        self.setLayout(layout)
        
        self.holdFigure = False

    def open_json_file(self):
        # Open file dialog to select a JSON file
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open JSON File", "", "JSON Files (*.json)")

        if file_path:
            self.open_json_file(file_path)
    
    def open_json_file(self, path):
        try:
            # Load the JSON file
            with open(path, 'r') as file:
                json_data = json.load(file)

            # Clear any previous errors
            self.error_label.clear()

            # Display the JSON in the tree view
            self.display_json(json_data)

        except (json.JSONDecodeError, FileNotFoundError) as e:
            # Show error message if there is an issue with the JSON
            self.error_label.setText(f"Error loading JSON: {e}")

    def display_json(self, json_data):
        # Create a model for the tree view
        model = QtGui.QStandardItemModel()
        model.setHorizontalHeaderLabels(['Key', 'Value'])

        # Start populating the tree view with JSON data
        root_item = model.invisibleRootItem()

        self.populate_tree(json_data, root_item)

        # Set the model on the tree view
        self.tree_view.setModel(model)

    def populate_tree(self, data, parent_item):
        if isinstance(data, dict):
            for key, value in data.items():
                key_item = QtGui.QStandardItem(key)
                if isprimitive(value):
                    parent_item.appendRow([key_item, QtGui.QStandardItem(str(value))])
                else:
                    parent_item.appendRow([key_item, QtGui.QStandardItem("")])
                    # Recursively add child nodes for the value
                    self.populate_tree(value, key_item)
        elif isinstance(data, list):
            for index, item in enumerate(data):
                index_item = QtGui.QStandardItem(f"[{index}]")
                if isprimitive(item):
                    parent_item.appendRow([index_item, QtGui.QStandardItem(str(item))])
                else:
                    parent_item.appendRow([index_item, QtGui.QStandardItem("")])
                    # Recursively add child nodes for the item
                    self.populate_tree(item, index_item)
        elif data is None:
            value_item = QtGui.QStandardItem("<None>")
            parent_item.setChild(0, value_item)
        else:
            # Base case: if it's a value, just display it
            value_item = QtGui.QStandardItem(str(data))
            parent_item.setChild(0, value_item)

    def plot_selected_as_monomial(self):
        self.plot_selected(power_eval)
    
    def plot_selected_as_chebyshev(self):
        self.plot_selected(cheb_eval)

    def plot_selected(self, eval_func):
        rootIdx = self.tree_view.selectedIndexes()[0]
        m = rootIdx.model()
        
        childCount = m.rowCount(rootIdx)
        
        plot_segments = []
        for child_i in range(childCount):
            childIdx = m.index(child_i, 0, rootIdx)
            coeffs_row = 0
            childRows = m.rowCount(childIdx)
            coeffs = []
            func_range = []
            func_domain = []
            for seg_i in range(childRows):
                segAttrIdx = m.index(seg_i, 0, childIdx)
                print("attribute: " + segAttrIdx.data())
                if "coeff" in segAttrIdx.data():
                    print("found " + segAttrIdx.data())
                    coeffCount = m.rowCount(segAttrIdx)
                    print("coeff count: " + str(coeffCount))
                    
                    for coeff_i in range(coeffCount):
                        coeffs.append(float(m.index(coeff_i, 1, segAttrIdx).data()))
                if "range" in segAttrIdx.data():
                    func_range.append(float(m.index(0,1,segAttrIdx).data()))
                    func_range.append(float(m.index(1,1,segAttrIdx).data()))
                if "domain" in segAttrIdx.data():
                    func_domain.append(float(m.index(0,1,segAttrIdx).data()))
                    func_domain.append(float(m.index(1,1,segAttrIdx).data()))
            print("range: " + str(func_range))
            print("domain: " + str(func_domain))
            print("coeffs: " + str(coeffs))

            points = max(10,int((func_domain[1]-func_domain[0])/0.05))
            xs = np.linspace(func_domain[0], func_domain[1], points)
            ys = eval_func(xs, coeffs, func_domain, func_range)
            plot_segments.append(go.Scatter(x=xs,y=ys,mode='lines'))
        
        # Create the layout
        layout = go.Layout(
        title="",
            xaxis=dict(title="x"),
            yaxis=dict(title="y"),
            showlegend=True
        )

        if self.holdFigure:
            print("adding trace")
            self.curfig.add_traces(plot_segments)
            self.curfig.show()
        else:
            print("creating new figure")
            self.curfig = go.Figure(data=plot_segments, layout=layout)
            self.curfig.show()
        # print("plotting: " + str(segments));