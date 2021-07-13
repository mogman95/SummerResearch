from PyQt5 import QtWidgets, uic
import sys
import os
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import sys


"""Analysis Window"""
class seqAnalysisWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(seqAnalysisWindow, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi('seq_analysis_window.ui', self) # Load the .ui file

        """ Connect functions """

        self.hide() # Show the GUI

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi('davids_main.ui', self) # Load the .ui file

        """ Connect functions """
        self.select_file_btn.clicked.connect(self.selectFileDialog)
        self.analyze_btn.clicked.connect(self.open_seq_window)
        self.AnalysisWindow = seqAnalysisWindow()
        self.AnalysisWindow.test_btn.clicked.connect(self.go_back)
#        self.AnalysisWindow.hide()
        self.show() # Show the GUI

    def selectFileDialog(self):
        file_filter = 'File Type(*);; FASTA (*.fasta *.fna);; Text (*.txt);; Data (*.xlsx *.xls *.csv *.dat)'
        response = QFileDialog.getOpenFileNames(
            # parent=self, ### This crashes the program if functioning
            caption='Select a data file',
            directory=os.getcwd(), ### gets "current working directory" (idk what that means)
            filter=file_filter
        )
        file_path = response[0][0] #response returns a list inside a duple => [0][0] indexes the file path out
        print(file_path)
        self.filepath_edit.setText(file_path)

    def open_seq_window(self):
        self.hide()
        self.AnalysisWindow.show()
        if self.strands_box.isChecked():
            self.AnalysisWindow.nucfreq_analysis_label.setHidden(True)
            self.AnalysisWindow.nucfreq_analysis_txtbrwsr.setHidden(True)

    def go_back(self):
        self.AnalysisWindow.hide()
        self.show()
            
app = QtWidgets.QApplication(sys.argv) # Create an instance of QtWidgets.QApplication
window = Ui() # Create an instance of our class
app.exec_() # Start the application