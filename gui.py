#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gui.py (PyQt5 版本)

提供一个基于 PyQt5 的图形界面，用于填写 core.py 所需的参数并启动整个流程：
  - 参考 FASTA 文件选择
  - 工作目录（输出）选择
  - 线程数输入
  - 样本测序文件夹选择
  - 点击“运行”后，后台调用 core.py，并把日志实时输出到界面下方的文本框

使用前提：
  1. 已安装 PyQt5：`pip install PyQt5`
  2. core.py 及所需的 GeneMiner2 原始脚本应与本 gui.py 放在同一目录
  3. 每个样本目录下需包含双端文件：<sample_name>_1.fastq 以及 <sample_name>_2.fastq （或同名 .gz 压缩）
"""

import sys
import os
import subprocess
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton, QListWidget,
    QTextEdit, QFileDialog, QHBoxLayout, QVBoxLayout, QMessageBox, QSpinBox
)
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtGui import QTextCursor, QColor

class WorkerThread(QThread):
    """
    后台线程用于执行 core.py，并实时把标准输出／错误输出发送回 GUI。
    """
    output_signal = pyqtSignal(str)   # 正常日志
    error_signal  = pyqtSignal(str)   # 错误日志
    finished_signal = pyqtSignal(int) # 退出码

    def __init__(self, cmd, workdir=None):
        super().__init__()
        self.cmd = cmd
        self.workdir = workdir

    def run(self):
        """
        启动子进程，循环读取 stdout/stderr，并把行内容通过信号发送给主线程。
        """
        try:
            # 使用 Popen 以实时读取输出
            process = subprocess.Popen(
                self.cmd,
                cwd=self.workdir,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                bufsize=1,
                universal_newlines=True,
            )
        except Exception as e:
            self.error_signal.emit(f"无法启动进程：{e}\n")
            self.finished_signal.emit(-1)
            return

        # 交替读取 stdout 和 stderr
        while True:
            out_line = process.stdout.readline()
            err_line = process.stderr.readline()

            if out_line:
                self.output_signal.emit(out_line)
            if err_line:
                self.error_signal.emit(err_line)

            if out_line == "" and err_line == "" and process.poll() is not None:
                break

        exit_code = process.wait()
        self.finished_signal.emit(exit_code)


class GeneMinerGUI(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("CPMiner2 GUI 主控 (PyQt5)")
        self.resize(800, 600)

        # 存储样本路径
        self.sample_dir = ""

        # 创建界面
        self._init_ui()

        # 用来保存 WorkerThread 对象的引用，防止被垃圾回收
        self.worker = None

    def _init_ui(self):
        """创建并排列所有控件"""

        # —— 参考 FASTA —— 
        lbl_ref = QLabel("参考 FASTA（multi-FASTA）：")
        self.edit_ref = QLineEdit()
        btn_ref_browse = QPushButton("浏览")
        btn_ref_browse.clicked.connect(self._browse_ref)

        hbox_ref = QHBoxLayout()
        hbox_ref.addWidget(lbl_ref)
        hbox_ref.addWidget(self.edit_ref)
        hbox_ref.addWidget(btn_ref_browse)

        # —— 工作目录 —— 
        lbl_outdir = QLabel("工作目录（输出）：")
        self.edit_outdir = QLineEdit()
        btn_outdir_browse = QPushButton("浏览")
        btn_outdir_browse.clicked.connect(self._browse_outdir)

        hbox_outdir = QHBoxLayout()
        hbox_outdir.addWidget(lbl_outdir)
        hbox_outdir.addWidget(self.edit_outdir)
        hbox_outdir.addWidget(btn_outdir_browse)

        # —— 线程数 —— 
        lbl_threads = QLabel("线程数：")
        self.spin_threads = QSpinBox()
        self.spin_threads.setRange(1, 64)
        self.spin_threads.setValue(1)

        hbox_threads = QHBoxLayout()
        hbox_threads.addWidget(lbl_threads)
        hbox_threads.addWidget(self.spin_threads)
        hbox_threads.addStretch()

        # —— 样本文件夹选择 —— 
        lbl_sample_dir = QLabel("样本文件夹：")
        self.edit_sample_dir = QLineEdit()
        btn_sample_browse = QPushButton("浏览")
        btn_sample_browse.clicked.connect(self._browse_sample_dir)

        hbox_sample_input = QHBoxLayout()
        hbox_sample_input.addWidget(lbl_sample_dir)
        hbox_sample_input.addWidget(self.edit_sample_dir)
        hbox_sample_input.addWidget(btn_sample_browse)

        # —— 运行 / 取消 按钮 —— 
        btn_run = QPushButton("运行")
        btn_run.clicked.connect(self._on_run)

        btn_quit = QPushButton("退出")
        btn_quit.clicked.connect(self.close)

        hbox_run = QHBoxLayout()
        hbox_run.addWidget(btn_run)
        hbox_run.addWidget(btn_quit)
        hbox_run.addStretch()

        # —— 日志输出区 —— 
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)

        # —— 整体布局 —— 
        main_layout = QVBoxLayout()
        main_layout.addLayout(hbox_ref)
        main_layout.addLayout(hbox_outdir)
        main_layout.addLayout(hbox_threads)
        main_layout.addSpacing(10)
        main_layout.addLayout(hbox_sample_input)
        main_layout.addLayout(hbox_run)
        main_layout.addWidget(self.txt_log)

        self.setLayout(main_layout)

    def _browse_ref(self):
        """弹出对话框选择参考 FASTA 文件"""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "选择参考 FASTA 文件",
            "",
            "FASTA Files (*.fasta *.fa *.fna *.gb *.genbank);;All Files (*)"
        )
        if path:
            self.edit_ref.setText(path)

    def _browse_outdir(self):
        """弹出对话框选择工作/输出目录"""
        path = QFileDialog.getExistingDirectory(
            self,
            "选择工作目录（若不存在会自动创建）",
            ""
        )
        if path:
            self.edit_outdir.setText(path)

    def _browse_sample_dir(self):
        """弹出对话框选择样本文件夹"""
        path = QFileDialog.getExistingDirectory(
            self,
            "选择包含样本的文件夹",
            ""
        )
        if path:
            self.edit_sample_dir.setText(path)

    def _on_run(self):
        """运行按钮被点击时，校验输入并调用 core.py"""
        ref     = self.edit_ref.text().strip()
        outdir  = self.edit_outdir.text().strip()
        threads = self.spin_threads.value()

        # 校验
        if not ref or not os.path.isfile(ref):
            QMessageBox.warning(self, "错误", "请选择合法存在的参考 FASTA 文件。")
            return
        if not outdir:
            QMessageBox.warning(self, "错误", "请选择工作目录。")
            return
        if not os.path.isdir(self.edit_sample_dir.text().strip()):
            QMessageBox.warning(self, "错误", "请选择合法的样本文件夹。")
            return

        # 检查并创建输出工作目录
        try:
            os.makedirs(outdir, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "错误", f"无法创建工作目录：{e}")
            return

        # 获取样本文件夹路径
        samples_arg = self.edit_sample_dir.text().strip()

        # 检查 core.py 是否存在
        script_dir = os.path.dirname(os.path.realpath(__file__))
        core_py    = os.path.join(script_dir, "core.py")
        if not os.path.isfile(core_py):
            QMessageBox.critical(self, "错误", f"找不到 core.py ({core_py})")
            return

        # 拼接命令行
        cmd = [
            sys.executable,
            core_py,
            "--ref", ref,
            "--samples", samples_arg,  # 传递文件夹路径
            "--outdir", outdir,
            "--threads", str(threads)
        ]

        # 禁用运行按钮，避免多次点击
        sender = self.sender()
        sender.setEnabled(False)

        # 清空日志区
        self.txt_log.clear()
        self.txt_log.append(f">>> 调用核心流程:\n  {' '.join(cmd)}\n")

        # 启动后台线程执行
        self.worker = WorkerThread(cmd, workdir=script_dir)
        self.worker.output_signal.connect(self._append_output)
        self.worker.error_signal.connect(self._append_error)
        self.worker.finished_signal.connect(self._proc_finished)
        self.worker.start()

    def _append_output(self, text):
        """把子进程 stdout 的内容追加到日志区"""
        self.txt_log.moveCursor(QTextCursor.End)
        self.txt_log.insertPlainText(text)

    def _append_error(self, text):
        """把子进程 stderr 的内容追加到日志区（以红色显示）"""
        self.txt_log.moveCursor(QTextCursor.End)
        self.txt_log.setTextColor(QColor("blue"))
        self.txt_log.insertPlainText(text)
        self.txt_log.setTextColor(QColor("black"))

    def _proc_finished(self, exit_code):
        """子进程结束后恢复按钮、输出提示"""
        QMessageBox.information(
            self, "提示",
            f"流程完成，退出码：{exit_code}\n请查看日志了解详细信息。"
        )
        # 重新启用“运行”按钮
        for w in self.findChildren(QPushButton):
            if w.text() == "运行":
                w.setEnabled(True)
                break


def main():
    app = QApplication(sys.argv)
    gui = GeneMinerGUI()
    gui.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
