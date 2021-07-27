"""Imports"""
import os
import sys
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets, uic
import analysisFunction as af # Imports a python file that does the analysis calculations

"""No File Dialog"""
class nofile_Dialog(QtWidgets.QDialog):        ### This dialog is presented when an invalid file/file path is entered
    def __init__(self):                        # Initialize class
        super(nofile_Dialog, self).__init__()  # Call the inherited classes __init__ method
        uic.loadUi('nofile_dialog.ui', self)   # Load the .ui file
        self.ok_btn.clicked.connect(self.hide) # Hide the dialog when the 'Ok' button is pressed

"""Help Dialog"""
class help_Dialog(QtWidgets.QDialog):          ### This dialog is presented when the help section in the analysis window is clicked
    def __init__(self):                        # Initialize class
        super(help_Dialog, self).__init__()    # Call the inherited classes __init__ method
        uic.loadUi('help_dialog.ui', self)     # Load the .ui file
        self.ok_btn.clicked.connect(self.hide) # Hide the dialog when the 'Ok' button is pressed

"""Analysis Window"""
class seqAnalysisWindow(QtWidgets.QMainWindow):    ### Window that presents the analysis results on the selected sequence
    def __init__(self):                            # Initialize class
        super(seqAnalysisWindow, self).__init__()  # Call the inherited classes __init__ method
        uic.loadUi('seq_analysis_window.ui', self) # Load the .ui file
        ### Connect to help dialog
        self.option_help = help_Dialog()                                      # Initialize help dialog
        self.actionFile_Options_Help.triggered.connect(self.option_help.show) # When the "File Options Help" tab is clicked, show the help dialog

    '''Initialize Functions'''
    ### Presents error message in the sequence text browser
    def nuc_show_error(self):
        error_msg = '!!! ERROR: Invalid Input !!!\nMake sure your input only includes positive integers no longer than the length of the sequence and is formated as:\n\"[starting nucleotide] to [ending nucleotide]\"'
        self.seq_analysis_txtbrwsr.setText(error_msg) # Set the text in the seq. text browser to the error msg.

    ### Sets the start_txt and end_txt variables from the respective edit lines
    def start_stop(self):
        ### Take the text in the edit line, strip it of spaces, save the result in the respective variable
        self.start_txt = self.show_nuc_start_edit.text().replace(' ','')
        self.end_txt = self.show_nuc_end_edit.text().replace(' ','')

        try:    # Attempt the following code
            ### Set the processed text to an integer object
            self.start = int(self.start_txt)
            self.end = int(self.end_txt)
            self.show_nucs_btn.setEnabled(True) # Enable the 'show' btn
        except: # If there is an error, run this code instead
            self.nuc_show_error()                # Present error msg
            self.show_nucs_btn.setDisabled(True) # Disable the 'show' btn

    ### Presents the section, or whole, of the sequence to the seq. text browser
    def seq_(self, sequence_):
        self.seq_analysis_txtbrwsr.setText(sequence_[self.start - 1 : self.end]) # Set the text browser text to what ever section of the sequence is specified by the start and end inputs

    ### Presents the section, or whole, of the complement to the seq. text browser
    def comp_(self, comp_seq_):
        self.seq_analysis_txtbrwsr.setText(comp_seq_[self.start - 1 : self.end])

    ### Presents the section, or whole, of the reverse complement to the seq. text browser
    def rev_comp_(self, rev_comp_, length_):
        self.seq_analysis_txtbrwsr.setText(rev_comp_[int(length_) - self.end : int(length_) - (self.start-1)])

    ### Checks which of the radio btns are selected and presents the appropriate sequence
    def checked_show(self, length_, sequence_, comp_seq_, rev_comp_):
        try:
            if (self.end - self.start) < 0 or self.end > int(length_): # If the start and end order is flipped OR the end is farther than the length of the sequence...
                raise Exception()                                      # Present error, skip code
            ### If a certain radio btn is selected, display its respective sequence
            if self.seq_radiobtn.isChecked():
                self.seq_(sequence_)
            elif self.comp_radiobtn.isChecked():
                self.comp_(comp_seq_)
            elif self.rev_comp_radiobtn.isChecked():
                self.rev_comp_(rev_comp_, length_)
        except:
            self.nuc_show_error()

    ### Resets the nucleotide selection edit lines and text browser
    def seq_reset(self, length_, sequence_, comp_seq_, rev_comp_):
        self.start = 1
        self.end = int(length_) # Make the length of the sequence an integer

        ### Set the text in the edit lines to the reset start and end values
        self.show_nuc_start_edit.setText(str(self.start))
        self.show_nuc_end_edit.setText(str(self.end))

        ### Avoids crashing program by overloading the text browser
        if int(length_) < 100000:                                       # If the length of the sequence is less than 100,000 nucleotides...
            self.checked_show(length_, sequence_, comp_seq_, rev_comp_) # Present appropriate sequence

    ### Hides all of the analysis widgets
    def hide_widgets(self):
        self.nucfreq_analysis_label.hide()
        self.nucfreq_analysis_txtbrwsr.hide()
        self.dinucfreq_analysis_label.hide()
        self.dinucfreq_analysis_txtbrwsr.hide()
        self.gc_cont_analysis_label.hide()
        self.gc_cont_analysis_txtbrwsr.hide()
        self.length_analysis_label.hide()
        self.length_analysis_txtbrwsr.hide()
        self.seq_view_analysis_label.hide()
        self.seq_radiobtn.hide()
        self.rev_comp_radiobtn.hide()
        self.comp_radiobtn.hide()
        self.whichseq_analysis_label.hide()
        self.show_nuc_label.hide()
        self.show_nuc_end_edit.hide()
        self.show_nuc_start_edit.hide()
        self.to_label.hide()
        self.seq_analysis_txtbrwsr.hide()
        self.reset_btn.hide()
        self.show_nucs_btn.hide()

"""Main Window"""
class analyseqMain(QtWidgets.QMainWindow):
    def __init__(self):                             # Initialize class
        super(analyseqMain, self).__init__()        # Call the inherited classes __init__ method
        uic.loadUi('analyseq-mainwindow.ui', self)  # Load the .ui file

        self.logo_label.setPixmap(QPixmap(u"Analyseq-logo.png")) # Sets the logo

        '''Connect functions'''
        self.select_file_btn.clicked.connect(self.selectFileDialog) # Prompts the user to select a file after the 'select file' button is pressed
        self.analyze_btn.clicked.connect(self.open_seq_window)      # Presents the analysis window after the 'analyze' button is pressed

        self.AnalysisWindow = seqAnalysisWindow()                   # Initializes the analysis window
        self.AnalysisWindow.hide_widgets()                          # Hides the analysis widgets
        self.AnalysisWindow.back_btn.clicked.connect(self.go_back)  # Presents the previous main window and hides the analysis window after the 'go back' button is presseed
        self.check_all_btn.clicked.connect(self.checkall)           # Checks all of the analysis options when the 'check all' button is pressed

        self.show() # Show the GUI

    '''Initialize Functions'''
    ### Prompts file selection
    def selectFileDialog(self):
        file_filter = 'FASTA (*.fasta *.fna)' # File filter has the label 'FASTA' and only calls .fasta and .fna files
        ### Prompts file selection and saves the file path and file type in a tuple
        response = QFileDialog.getOpenFileNames(caption='Select a data file', # Sets the title of the window
                                                directory=os.getcwd(),        # Sets the starting directory as the working directory
                                                filter=file_filter)           # Sets the file filter
        try:
            self.file_path = response[0][0]            # response returns a list inside a tuple => [0][0] indexes the file path out
            self.filepath_edit.setText(self.file_path) # Sets the file path edit line text to the file path
        except:
            self.filepath_edit.setText('Please select a file') # If no file is selected the file path edit line displays 'Please select a file'

    ### Presents the analysis window with all of the selected analyses
    def open_seq_window(self):
        try:
            if os.path.exists(self.file_path): # If the the selected file path is valid...
                self.hide()                    # Hide the main window
                self.AnalysisWindow.show()     # Show the analysis window

                data = af.main(self.file_path) # Analyzes the selected file and saves the data as a list in a variable
                ### Index each analysis out
                sequence = data[0]
                comp_seq = data[1]
                rev_comp = data[2]
                GC_cont = data[3]
                length = data[4]
                nuc_results = data[5]
                dinuc_results = data[6]

                self.AnalysisWindow.seq_reset(length, sequence, comp_seq, rev_comp) # Primes the sequence presentation analysis

                ### When the text in the start and end edit lines is changed, process the start and end variables
                self.AnalysisWindow.show_nuc_start_edit.textChanged.connect(lambda *args: self.AnalysisWindow.start_stop())
                self.AnalysisWindow.show_nuc_end_edit.textChanged.connect(lambda *args: self.AnalysisWindow.start_stop())

                self.AnalysisWindow.show_nucs_btn.clicked.connect(lambda *args: self.AnalysisWindow.checked_show(length, sequence, comp_seq, rev_comp)) # Show the desired section of the sequence when the 'show' btn is pressed
                self.AnalysisWindow.reset_btn.clicked.connect(lambda *args: self.AnalysisWindow.seq_reset(length, sequence, comp_seq, rev_comp))        # Reset sequence shown when the 'reset' btn is pressed

                self.AnalysisWindow.seq_radiobtn.toggled.connect(lambda *args: self.AnalysisWindow.seq_(sequence))           # Show the sequence when the sequence radio btn is selected
                self.AnalysisWindow.comp_radiobtn.toggled.connect(lambda *args: self.AnalysisWindow.comp_(comp_seq))         # Show the complement when the complement radio btn is selected
                self.AnalysisWindow.rev_comp_radiobtn.toggled.connect(lambda *args: self.AnalysisWindow.rev_comp_(rev_comp, length)) # Show the reverse complement when the reverse complement radio btn is selected

                self.AnalysisWindow.actionNew.triggered.connect(lambda *args: self.new_analysis()) # Presents a blank main window when the 'New analysis' tab is pressed
                self.AnalysisWindow.actionsave.triggered.connect(lambda *args: self.Save_(length, GC_cont, nuc_results, dinuc_results, sequence, comp_seq, rev_comp))      # Saves a txt file of all of the analyses when the 'Save' tab is pressed
                self.AnalysisWindow.actionSave_As.triggered.connect(lambda *args: self.save_as(length, GC_cont, nuc_results, dinuc_results, sequence, comp_seq, rev_comp)) # Prompt the user to save a custom file name and type of the selected analyses
                self.AnalysisWindow.actionClose.triggered.connect(lambda *args: self.go_back()) # Same action as the 'Go Back' btn when the 'Close' tab is pressed
                self.AnalysisWindow.actionQuit.triggered.connect(lambda *args: sys.exit())      # Closes the application when the 'Quit' tab is pressed

                ### If an analysis is checked in the main window, show its respective analysis in the analysis window
                if self.nuc_freq_box.isChecked():
                    self.AnalysisWindow.nucfreq_analysis_txtbrwsr.setText(nuc_results)
                    self.AnalysisWindow.nucfreq_analysis_label.show()
                    self.AnalysisWindow.nucfreq_analysis_txtbrwsr.show()
                if self.dinuc_freq_box.isChecked():
                    self.AnalysisWindow.dinucfreq_analysis_txtbrwsr.setText(dinuc_results)
                    self.AnalysisWindow.dinucfreq_analysis_label.show()
                    self.AnalysisWindow.dinucfreq_analysis_txtbrwsr.show()
                if self.gc_cont_box.isChecked():
                    self.AnalysisWindow.gc_cont_analysis_txtbrwsr.setText(GC_cont + '%')
                    self.AnalysisWindow.gc_cont_analysis_label.show()
                    self.AnalysisWindow.gc_cont_analysis_txtbrwsr.show()
                if self.length_box.isChecked():
                    self.AnalysisWindow.length_analysis_txtbrwsr.setText(length)
                    self.AnalysisWindow.length_analysis_label.show()
                    self.AnalysisWindow.length_analysis_txtbrwsr.show()
                if self.strands_box.isChecked():
                    self.AnalysisWindow.seq_view_analysis_label.show()
                    self.AnalysisWindow.seq_radiobtn.show()
                    self.AnalysisWindow.rev_comp_radiobtn.show()
                    self.AnalysisWindow.comp_radiobtn.show()
                    self.AnalysisWindow.whichseq_analysis_label.show()
                    self.AnalysisWindow.show_nuc_label.show()
                    self.AnalysisWindow.show_nuc_end_edit.show()
                    self.AnalysisWindow.show_nuc_start_edit.show()
                    self.AnalysisWindow.to_label.show()
                    self.AnalysisWindow.reset_btn.show()
                    self.AnalysisWindow.show_nucs_btn.show()
                    ### The following code prevents the application from crashing from overloading the text browser with a sequence that is too long
                    if self.AnalysisWindow.seq_radiobtn.isChecked() and int(length) < 100000: # If the sequence radio btn is selected (which it is by default), and the length of the sequence is less than 100,000...
                        self.AnalysisWindow.seq_(sequence)                                    # Set the sequence in the text browser

                    if int(length) > 100000: # If the length of the sequence is greater than 100,000...
                        ### Presents a message informing the user that the sequence is long
                        self.AnalysisWindow.seq_analysis_txtbrwsr.setText(f'Sequence is {length} nucleotides long.\nIt is recommended to show less than 100,000 nucleotides at a time')
                        self.AnalysisWindow.reset_btn.setDisabled(True)         # Disable the 'reset' btn
                        self.AnalysisWindow.show_nuc_end_edit.setText('100000') # Present only the first 100,000 nucleotides by default
                        self.AnalysisWindow.seq_analysis_txtbrwsr.show()        # Show the text browser
                    else:
                        self.AnalysisWindow.seq_analysis_txtbrwsr.show() # Show the full sequence
            else:
                self.no_file_dialog() # Show the no file dialog telling the user to select a valid file
        except:
            self.no_file_dialog()

    ### Initializes and presents the no file dialog
    def no_file_dialog(self):
        self.dialog = nofile_Dialog()
        self.dialog.show()

    ### Hides the analysis window and presents the previous main window
    def go_back(self):
        self.AnalysisWindow.hide()
        self.AnalysisWindow.hide_widgets() # Resets the widgets incase the user selects different analyses
        self.show()

    ### Resets the main window, hides the analysis window, and presents the blank main window
    def new_analysis(self):
        self.filepath_edit.setText('')       # Clears the file path edit line
        self.check_all_btn.setChecked(False) # Unchecks the 'check all' btn if it is checked
        self.checkall()                      # Uncheck all of the analyses
        self.AnalysisWindow.hide()
        self.AnalysisWindow.hide_widgets()
        self.show()

    ### Checks or unchecks all of the analyses based on the checked state of the 'check all' btn
    def checkall(self):
        if self.check_all_btn.isChecked():
            self.nuc_freq_box.setCheckState(Qt.Checked)
            self.dinuc_freq_box.setCheckState(Qt.Checked)
            self.gc_cont_box.setCheckState(Qt.Checked)
            self.length_box.setCheckState(Qt.Checked)
            self.strands_box.setCheckState(Qt.Checked)
        else:
            self.nuc_freq_box.setCheckState(Qt.Unchecked)
            self.dinuc_freq_box.setCheckState(Qt.Unchecked)
            self.gc_cont_box.setCheckState(Qt.Unchecked)
            self.length_box.setCheckState(Qt.Unchecked)
            self.strands_box.setCheckState(Qt.Unchecked)

    ### Saves a txt file of all of the analyses
    def Save_(self, length_, GC_cont_, nuc_results_, dinuc_results_, sequence_, comp_seq_, rev_comp_):
        output_file = "QuickSaveAnalysis.txt"
        with open(output_file, 'w+') as f:  # Open a new file, or override and old file, that can be written to as 'f'
            ### Writes each of the following lines to the text file in order
            f.write("\nLength:\n" + length_)
            f.write("\nGC Content:\n" + GC_cont_)
            f.write("\nNucleotide Frequencies:\n" + nuc_results_)
            f.write("\nDinucleotide Frequencies:\n" + dinuc_results_)
            f.write("\nSequence:\n" + sequence_)
            f.write("\n\nComplement:\n" + comp_seq_)
            f.write("\n\nReverse Complement:\n" + rev_comp_)

    ### Prompt the user to save a custom file name and type of the selected analyses
    def save_as(self, length_, GC_cont_, nuc_results_, dinuc_results_, sequence_, comp_seq_, rev_comp_):
        try:
            file_filter = 'TEXT (*.txt);; CSV (*.csv);; HTLM (*.htlm)' # The different file types the user has to choose from
            ### Prompts the 'save file' window and saves the file path and type in a tuple
            name = QFileDialog.getSaveFileName(self, caption='Save As...',  # Sets the window title
                                                     directory=os.getcwd(), # Sets the starting directory as the working directory
                                                     filter=file_filter     # Set the file types to choose from
                                               )[0]                         # Index the file path out so that the file path is saved in 'name' as a string
            with open(name, 'w') as f:
                ### Checks if the analysis is checked and, if it is, adds it to the file
                if self.nuc_freq_box.isChecked():
                    f.write("\nNucleotide Frequencies:\n" + nuc_results_)
                if self.dinuc_freq_box.isChecked():
                    f.write("\nDinucleotide Frequencies:\n" + dinuc_results_)
                if self.gc_cont_box.isChecked():
                    f.write("\nGC Content:\n" + GC_cont_)
                if self.length_box.isChecked():
                    f.write("\nLength:\n" + length_)
                ### All 3 strands are included for simplicity
                if self.strands_box.isChecked():
                    f.write("\n\nSequence:\n" + sequence_)
                    f.write("\n\nComplement:\n" + comp_seq_)
                    f.write("\n\nReverse Complement:\n" + rev_comp_)
        except:       # If no file is made (an exception is raised)...
            'Go back' # Do nothing and go back

app = QtWidgets.QApplication(sys.argv)  # Create an instance of QtWidgets.QApplication
window = analyseqMain()                 # Create an instance of our class
app.exec_()                             # Start the application