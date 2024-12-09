import numpy as np
import pandas as pd
import os, sys
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import gridspec
from MATS import linelistdata
parent_dir = os.path.abspath('../../')
sys.path.append(parent_dir)
import MATS
import seaborn as sns
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("poster")
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QSpinBox, QLabel, QPushButton
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon
from MATS.hapi import *
# 模拟的光谱数据生成函数
def simulate_spectrum(temperature, pressure, line_shape):
    # 假设这里是简单的模拟方法，实际可以替换为更复杂的模型
    # 模拟的波数范围，单位为 cm^-1
    wavenumbers = np.linspace(2000, 2200, 500)  # 2000-2200 cm^-1

    # 基本的高斯线型模拟
    if line_shape == 'Gaussian':
        # 基础参数
        gamma = 0.1  # 宽度
        intensity = np.exp(-0.5 * ((wavenumbers - 2100) / gamma) ** 2)  # 高斯分布的强度

    # 简单的压力和温度对模拟的影响
    intensity *= (pressure / 1.0) * (temperature / 300.0)  # 假设的比例关系

    return wavenumbers, intensity

# 主窗口类
class SpectrumSimulator(QWidget):
    def __init__(self):
        super().__init__()

        # 设置窗口
        self.setWindowTitle('光谱模拟软件')
        self.setGeometry(100, 100, 600, 400)
        self.setWindowIcon(QIcon("icon.png"))  # 你可以替换为你自己的图标
        db_begin('data')
        #fetch('CO2', 2, 1, 2000, 2100)

        # 初始化UI
        self.initUI()

    def initUI(self):
        # 主布局

        main_layout = QVBoxLayout()
        # 创建温度输入框
        temp_label = QLabel('温度 (K):')
        self.temp_spinbox = QSpinBox()
        self.temp_spinbox.setRange(150, 5000)  # 温度范围
        self.temp_spinbox.setValue(300)  # 默认值为300K

        # 创建压力输入框
        pressure_label = QLabel('压力 (atm):')
        self.pressure_spinbox = QSpinBox()
        self.pressure_spinbox.setRange(0, 1000)  # 压力范围
        self.pressure_spinbox.setValue(100)  # 默认值为1 atm

        # 创建线型选择框
        line_shape_label = QLabel('选择线型:')
        self.line_shape_combobox = QComboBox()
        self.line_shape_combobox.addItems(['Gaussian', 'Lorentzian','VP','SDVP','HT'])

        # 创建模拟按钮
        self.simulate_button = QPushButton('生成光谱')
        self.simulate_button.clicked.connect(self.generate_spectrum)

        # 设置布局
        form_layout = QVBoxLayout()
        form_layout.addWidget(temp_label)
        form_layout.addWidget(self.temp_spinbox)
        form_layout.addWidget(pressure_label)
        form_layout.addWidget(self.pressure_spinbox)
        form_layout.addWidget(line_shape_label)
        form_layout.addWidget(self.line_shape_combobox)
        form_layout.addWidget(self.simulate_button)

        main_layout.addLayout(form_layout)

        # 图表显示
        self.fig, self.ax = plt.subplots()
        self.canvas = self.fig.canvas

        # 将图表放入布局中
        self.canvas_widget = QWidget(self)
        self.canvas_layout = QVBoxLayout(self.canvas_widget)
        self.canvas_layout.addWidget(self.canvas)
        main_layout.addWidget(self.canvas_widget)
        # toolbar = NavigationToolbar(self.canvas)
        # main_layout.addWidget(toolbar)
        # 设置窗体的主布局
        self.setLayout(main_layout)

    # 线型选择函数
    def line_shape_changed(self,tablename):
        if self.line_shape == 'Gaussian':
            self.wavenumbers, self.intensity = absorptionCoefficient_Doppler(SourceTables=tablename, Diluent={'air': 1.0})
        elif self.line_shape == 'Lorentzian':
            self.wavenumbers, self.intensity = absorptionCoefficient_Lorentz(SourceTables=tablename, Diluent={'air': 1.0})
        elif self.line_shape == 'VP':
            self.wavenumbers, self.intensity = absorptionCoefficient_Voigt(SourceTables=tablename, Diluent={'air': 1.0})
        elif self.line_shape == 'SDVP':
            self.wavenumbers, self.intensity = absorptionCoefficient_SDVoigt(SourceTables=tablename, Diluent={'air': 1.0})
        elif self.line_shape == 'HT':
            self.wavenumbers, self.intensity = absorptionCoefficient_HT(SourceTables=tablename, Diluent={'air': 1.0})
        else:
            raise ValueError('Invalid line shape')
    # 生成光谱并绘制
    def generate_spectrum(self):
        # 获取输入的温度、压力和线型
        self.temperature = self.temp_spinbox.value()
        self.pressure = self.pressure_spinbox.value()
        self.line_shape = self.line_shape_combobox.currentText()

        # 使用模拟函数生成光谱
        tablename = 'CO2'
        #wavenumbers, intensity = absorptionCoefficient_Lorentz(SourceTables='CO2', Diluent={'air': 1.0})
        #wavenumbers, intensity = simulate_spectrum(temperature, pressure, line_shape)
        self.line_shape_changed(tablename)
        # 清空当前图表
        self.ax.clear()

        # 绘制新的光谱图
        self.ax.plot(self.wavenumbers, self.intensity, label=f'{self.line_shape} Line', color='blue')
        self.ax.set_xlabel('波数 (cm^-1)',fontdict={'fontname': 'SimHei'})
        self.ax.set_ylabel('强度',fontdict={'fontname': 'SimHei'})
        self.ax.set_title(f'模拟光谱 (T={self.temperature}K, P={self.pressure}atm)',fontdict={'fontname': 'SimHei'})
        self.ax.legend()

        # 更新图表
        #self.fig.tight_layout()
        self.canvas.draw()

# 主程序入口
def main():
    app = QApplication(sys.argv)
    window = SpectrumSimulator()
    window.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()