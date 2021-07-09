import sys
import os
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets



"""Startup Window"""
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(249*3, 292*3)
        MainWindow.setMinimumSize(QSize(249*3, 292*3))
        MainWindow.setMaximumSize(QSize(249*3, 292*3))

        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")

        self.nuc_freq_box = QCheckBox(self.centralwidget)
        self.nuc_freq_box.setObjectName(u"nuc_freq_box")
        self.nuc_freq_box.setGeometry(QRect(9*3, 124*3, 228*3, 17*3))

        self.dinuc_freq_box = QCheckBox(self.centralwidget)
        self.dinuc_freq_box.setObjectName(u"dinuc_freq_box")
        self.dinuc_freq_box.setGeometry(QRect(9*3, 147*3, 142*3, 18*3))

        self.gc_cont_box = QCheckBox(self.centralwidget)
        self.gc_cont_box.setObjectName(u"gc_cont_box")
        self.gc_cont_box.setGeometry(QRect(9*3, 171*3, 228*3, 17*3))

        self.length_box = QCheckBox(self.centralwidget)
        self.length_box.setObjectName(u"length_box")
        self.length_box.setGeometry(QRect(9*3, 194*3, 55*3, 17*3))

        self.check_all_btn = QToolButton(self.centralwidget)
        self.check_all_btn.setCheckable(True)
        self.check_all_btn.setObjectName(u"check_all_btn")
        self.check_all_btn.setGeometry(QRect(9*3, 101*3, 58*3, 17*3))
        self.check_all_btn.setCursor(QCursor(Qt.ArrowCursor))
        self.check_all_btn.clicked.connect(self.checkall)

        self.strands_box = QCheckBox(self.centralwidget)
        self.strands_box.setObjectName(u"strands_box")
        self.strands_box.setGeometry(QRect(9*3, 217*3, 228*3, 17*3))

        self.analyseq_logo = QLabel(self.centralwidget)
        self.analyseq_logo.setObjectName(u"analyseq_logo")
        self.analyseq_logo.setGeometry(QRect(-1*3, -1*3, 251*3, 41*3))
        self.analyseq_logo.setStyleSheet(u"background-image: url(:/newPrefix/Analyseq-logo.png);")
        self.analyseq_logo.setPixmap(QPixmap(u"C:/Users/Mogma/Downloads/Trinh Lab/SummerResearch/fasta Analysis/Analyseq-logo.png"))
        self.analyseq_logo.setScaledContents(True)

        self.select_file_btn = QToolButton(self.centralwidget)
        self.select_file_btn.clicked.connect(self.selectFileDialog)
        self.select_file_btn.setObjectName(u"select_file_btn")
        self.select_file_btn.setGeometry(QRect(10*3, 60*3, 62*3, 21*3))

        self.select_analysis_label = QLabel(self.centralwidget)
        self.select_analysis_label.setObjectName(u"select_analysis_label")
        self.select_analysis_label.setGeometry(QRect(9*3, 81*3, 222*3, 16*3))

        self.select_file_label = QLabel(self.centralwidget)
        self.select_file_label.setObjectName(u"select_file_label")
        self.select_file_label.setGeometry(QRect(9*3, 40*3, 161*3, 16*3))

        self.analysisWindow = seqAnalysisWindow()
        self.analyze_btn = QToolButton(self.centralwidget)
        self.analyze_btn.setObjectName(u"analyze_btn")
        self.analyze_btn.setGeometry(QRect(10*3, 240*3, 231*3, 31*3))
        # self.analyze_btn.clicked.connect(self.analysisWindow)

        self.filepath_edit = QLineEdit(self.centralwidget)
        self.filepath_edit.setObjectName(u"filepath_edit")
        self.filepath_edit.setGeometry(QRect(80*3, 60*3, 161*3, 20*3))
        # self.filepath_edit.setText()

        MainWindow.setCentralWidget(self.centralwidget)

        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"Analyseq Main Window", None))
        self.nuc_freq_box.setText(QCoreApplication.translate("MainWindow", u"Nucleotide Frequencies", None))
        self.dinuc_freq_box.setText(QCoreApplication.translate("MainWindow", u"Dinucleotide Frequencies", None))
        self.gc_cont_box.setText(QCoreApplication.translate("MainWindow", u"GC Content", None))
        self.length_box.setText(QCoreApplication.translate("MainWindow", u"Length", None))
        self.check_all_btn.setText(QCoreApplication.translate("MainWindow", u"Check All", None))
        self.strands_box.setText(QCoreApplication.translate("MainWindow", u"Show Sequence Strands", None))
        self.select_file_btn.setText(QCoreApplication.translate("MainWindow", u"Select File", None))
        self.select_analysis_label.setText(QCoreApplication.translate("MainWindow", u"Select the analysis you would like to perform:", None))
        self.select_file_label.setText(QCoreApplication.translate("MainWindow", u"Select a fasta file to be analyzed:", None))
        self.analyze_btn.setText(QCoreApplication.translate("MainWindow", u"Analyze", None))
    # retranslateUi

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

    # def showWindow(self):
    #     analysisWindow = seqAnalysisWindow()
    #     self.analysisWindow.show()


"""Analysis Window"""
class seqAnalysisWindow(object):
    def setupUi(self, seqAnalysisWindow):
        if not seqAnalysisWindow.objectName():
            seqAnalysisWindow.setObjectName(u"seqAnalysisWindow")
        seqAnalysisWindow.resize(381, 568)
        self.actionOpen = QAction(seqAnalysisWindow)
        self.actionOpen.setObjectName(u"actionOpen")
        self.actionsave = QAction(seqAnalysisWindow)
        self.actionsave.setObjectName(u"actionsave")
        self.actionSave_As = QAction(seqAnalysisWindow)
        self.actionSave_As.setObjectName(u"actionSave_As")
        self.actionClose = QAction(seqAnalysisWindow)
        self.actionClose.setObjectName(u"actionClose")
        self.actionQuit = QAction(seqAnalysisWindow)
        self.actionQuit.setObjectName(u"actionQuit")
        self.actionNew = QAction(seqAnalysisWindow)
        self.actionNew.setObjectName(u"actionNew")
        self.actionClose_All_Analysis = QAction(seqAnalysisWindow)
        self.actionClose_All_Analysis.setObjectName(u"actionClose_All_Analysis")
        self.centralwidget = QWidget(seqAnalysisWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.formLayout = QFormLayout(self.centralwidget)
        self.formLayout.setObjectName(u"formLayout")


        self.seq_view_analysis_label = QLabel(self.centralwidget)
        self.seq_view_analysis_label.setObjectName(u"seq_view_analysis_label")

        self.formLayout.setWidget(8, QFormLayout.SpanningRole, self.seq_view_analysis_label)

        self.radiobtn_layout = QHBoxLayout()
        self.radiobtn_layout.setObjectName(u"radiobtn_layout")
        self.seq_radiobtn = QRadioButton(self.centralwidget)
        self.seq_radiobtn.setObjectName(u"seq_radiobtn")

        self.radiobtn_layout.addWidget(self.seq_radiobtn)

        self.comp_radiobtn = QRadioButton(self.centralwidget)
        self.comp_radiobtn.setObjectName(u"comp_radiobtn")

        self.radiobtn_layout.addWidget(self.comp_radiobtn)

        self.rev_comp_radiobtn = QRadioButton(self.centralwidget)
        self.rev_comp_radiobtn.setObjectName(u"rev_comp_radiobtn")

        self.radiobtn_layout.addWidget(self.rev_comp_radiobtn)


        self.formLayout.setLayout(9, QFormLayout.LabelRole, self.radiobtn_layout)

        self.whichseq_analysis_label = QLabel(self.centralwidget)
        self.whichseq_analysis_label.setObjectName(u"whichseq_analysis_label")

        self.formLayout.setWidget(10, QFormLayout.SpanningRole, self.whichseq_analysis_label)

        self.show_nuc_layout = QHBoxLayout()
        self.show_nuc_layout.setObjectName(u"show_nuc_layout")
        self.show_nuc_label = QLabel(self.centralwidget)
        self.show_nuc_label.setObjectName(u"show_nuc_label")

        self.show_nuc_layout.addWidget(self.show_nuc_label)

        self.show_nuc_start_edit = QLineEdit(self.centralwidget)
        self.show_nuc_start_edit.setObjectName(u"show_nuc_start_edit")

        self.show_nuc_layout.addWidget(self.show_nuc_start_edit)

        self.to_label = QLabel(self.centralwidget)
        self.to_label.setObjectName(u"to_label")

        self.show_nuc_layout.addWidget(self.to_label)

        self.show_nuc_end_edit = QLineEdit(self.centralwidget)
        self.show_nuc_end_edit.setObjectName(u"show_nuc_end_edit")

        self.show_nuc_layout.addWidget(self.show_nuc_end_edit)


        self.formLayout.setLayout(11, QFormLayout.LabelRole, self.show_nuc_layout)

        self.seq_analysis_txtbrwsr = QTextBrowser(self.centralwidget)
        self.seq_analysis_txtbrwsr.setObjectName(u"seq_analysis_txtbrwsr")

        self.formLayout.setWidget(12, QFormLayout.SpanningRole, self.seq_analysis_txtbrwsr)




        self.nucfreq_analysis_txtbrwsr = QTextBrowser(self.centralwidget)
        self.nucfreq_analysis_txtbrwsr.setObjectName(u"nucfreq_analysis_txtbrwsr")

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.nucfreq_analysis_txtbrwsr)

        self.dinucfreq_analysis_txtbrwsr = QTextBrowser(self.centralwidget)
        self.dinucfreq_analysis_txtbrwsr.setObjectName(u"dinucfreq_analysis_txtbrwsr")

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.dinucfreq_analysis_txtbrwsr)

        self.gc_cont_analysis_txtbrwsr = QTextBrowser(self.centralwidget)
        self.gc_cont_analysis_txtbrwsr.setObjectName(u"gc_cont_analysis_txtbrwsr")

        self.formLayout.setWidget(5, QFormLayout.LabelRole, self.gc_cont_analysis_txtbrwsr)

        self.length_analysis_txtbrwsr = QTextBrowser(self.centralwidget)
        self.length_analysis_txtbrwsr.setObjectName(u"length_analysis_txtbrwsr")

        self.formLayout.setWidget(7, QFormLayout.LabelRole, self.length_analysis_txtbrwsr)

        self.length_analysis_label = QLabel(self.centralwidget)
        self.length_analysis_label.setObjectName(u"length_analysis_label")

        self.formLayout.setWidget(6, QFormLayout.LabelRole, self.length_analysis_label)

        self.gc_cont_analysis_label = QLabel(self.centralwidget)
        self.gc_cont_analysis_label.setObjectName(u"gc_cont_analysis_label")

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.gc_cont_analysis_label)

        self.dinucfreq_analysis_label = QLabel(self.centralwidget)
        self.dinucfreq_analysis_label.setObjectName(u"dinucfreq_analysis_label")

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.dinucfreq_analysis_label)

        self.nucfreq_analysis_label = QLabel(self.centralwidget)
        self.nucfreq_analysis_label.setObjectName(u"nucfreq_analysis_label")

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.nucfreq_analysis_label)

        seqAnalysisWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(seqAnalysisWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 381, 22))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        seqAnalysisWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(seqAnalysisWindow)
        self.statusbar.setObjectName(u"statusbar")
        seqAnalysisWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuFile.addAction(self.actionNew)
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionsave)
        self.menuFile.addAction(self.actionSave_As)
        self.menuFile.addAction(self.actionClose)
        self.menuFile.addAction(self.actionClose_All_Analysis)
        self.menuFile.addAction(self.actionQuit)

        self.retranslateUi(seqAnalysisWindow)

        QMetaObject.connectSlotsByName(seqAnalysisWindow)
    # setupUi

    def retranslateUi(self, seqAnalysisWindow):
        seqAnalysisWindow.setWindowTitle(QCoreApplication.translate("seqAnalysisWindow", u"Sequence Analysis", None))
        self.actionOpen.setText(QCoreApplication.translate("seqAnalysisWindow", u"Open", None))
        self.actionsave.setText(QCoreApplication.translate("seqAnalysisWindow", u"Save", None))
        self.actionSave_As.setText(QCoreApplication.translate("seqAnalysisWindow", u"Save As...", None))
        self.actionClose.setText(QCoreApplication.translate("seqAnalysisWindow", u"Close Current Analysis", None))
        self.actionQuit.setText(QCoreApplication.translate("seqAnalysisWindow", u"Quit", None))
        self.actionNew.setText(QCoreApplication.translate("seqAnalysisWindow", u"New Analysis", None))
        self.actionClose_All_Analysis.setText(QCoreApplication.translate("seqAnalysisWindow", u"Close All Analysis", None))
        self.seq_view_analysis_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"Select nucleotide sequence view:", None))
        self.seq_radiobtn.setText(QCoreApplication.translate("seqAnalysisWindow", u"Main Sequence", None))
        self.comp_radiobtn.setText(QCoreApplication.translate("seqAnalysisWindow", u"Complement", None))
        self.rev_comp_radiobtn.setText(QCoreApplication.translate("seqAnalysisWindow", u"Reverse Complement", None))
        self.whichseq_analysis_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"TextLabel", None))
        self.show_nuc_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"Show nucleotides", None))
        self.to_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"to", None))
        self.length_analysis_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"Length", None))
        self.gc_cont_analysis_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"GC Content", None))
        self.dinucfreq_analysis_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"Dinucleotide Frequencies", None))
        self.nucfreq_analysis_label.setText(QCoreApplication.translate("seqAnalysisWindow", u"Nucleotide Frequencies", None))
        self.menuFile.setTitle(QCoreApplication.translate("seqAnalysisWindow", u"File", None))
        self.menuHelp.setTitle(QCoreApplication.translate("seqAnalysisWindow", u"Help", None))
    # retranslateUi




app = QApplication(sys.argv)
MainWindow1 = QMainWindow()
MainWindow2 = QMainWindow()
ui = Ui_MainWindow()
ui.setupUi(MainWindow1)
ui2 = seqAnalysisWindow()
ui2.setupUi(MainWindow2)
MainWindow1.show()
MainWindow2.show()
sys.exit(app.exec_())