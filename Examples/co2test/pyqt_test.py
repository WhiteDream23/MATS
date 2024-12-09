import sys
from PyQt6.QtWidgets import QApplication, QMainWindow
from PyQt6.uic.Compiler.qtproxies import QtWidgets

from untitled2 import Ui_MainWindow


class MyMainForm(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(MyMainForm, self).__init__()
        self.setupUi(self)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QApplication(sys.argv)
    myw = MyMainForm()
    myw.show()
    sys.exit(app.exec())