# [Part 1]：主程序入口 + 通用后台运行逻辑 + 主界面布局框架

import sys
import os
import subprocess
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog,
    QHBoxLayout, QVBoxLayout, QMessageBox, QSpinBox, QTabWidget, QCheckBox,
    QDoubleSpinBox, QSplitter, QSizePolicy
)
from PyQt5.QtCore import QThread, pyqtSignal, Qt
from PyQt5.QtGui import QTextCursor, QColor

class WorkerThread(QThread):
    output_signal = pyqtSignal(str)
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(int)

    def __init__(self, cmd, workdir=None):
        super().__init__()
        self.cmd = cmd
        self.workdir = workdir

    def run(self):
        try:
            process = subprocess.Popen(
                self.cmd, cwd=self.workdir, shell=False,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                bufsize=1, universal_newlines=True
            )
        except Exception as e:
            self.error_signal.emit(f"无法启动进程：{e}\n")
            self.finished_signal.emit(-1)
            return

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

class Function1Tab(QWidget):
    def __init__(self):
        super().__init__()

        self.worker = None
        self._init_ui()

    def _init_ui(self):
        # ===== 左边输入参数区域 =====
        param_layout = QVBoxLayout()

        lbl_ref = QLabel("参考 GenBank 文件：")
        self.edit_ref = QLineEdit()
        btn_ref = QPushButton("浏览")
        btn_ref.clicked.connect(self._browse_ref)
        hbox_ref = QHBoxLayout()
        hbox_ref.addWidget(self.edit_ref)
        hbox_ref.addWidget(btn_ref)
        param_layout.addWidget(lbl_ref)
        param_layout.addLayout(hbox_ref)

        lbl_outdir = QLabel("输出目录：")
        self.edit_outdir = QLineEdit()
        btn_outdir = QPushButton("浏览")
        btn_outdir.clicked.connect(self._browse_outdir)
        hbox_out = QHBoxLayout()
        hbox_out.addWidget(self.edit_outdir)
        hbox_out.addWidget(btn_outdir)
        param_layout.addWidget(lbl_outdir)
        param_layout.addLayout(hbox_out)

        lbl_threads = QLabel("线程数：")
        self.spin_threads = QSpinBox()
        self.spin_threads.setRange(1, 64)
        self.spin_threads.setValue(4)
        param_layout.addWidget(lbl_threads)
        param_layout.addWidget(self.spin_threads)

        lbl_samples = QLabel("样本目录：")
        self.edit_sample_dir = QLineEdit()
        btn_samples = QPushButton("浏览")
        btn_samples.clicked.connect(self._browse_sample_dir)
        hbox_samples = QHBoxLayout()
        hbox_samples.addWidget(self.edit_sample_dir)
        hbox_samples.addWidget(btn_samples)
        param_layout.addWidget(lbl_samples)
        param_layout.addLayout(hbox_samples)

        self.btn_run = QPushButton("运行")
        self.btn_run.clicked.connect(self._on_run)
        param_layout.addWidget(self.btn_run)
        param_layout.addStretch()

        # ===== 右边日志区域 =====
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)

        # ===== 整体左右布局 =====
        splitter = QSplitter(Qt.Horizontal)
        param_widget = QWidget()
        param_widget.setLayout(param_layout)
        splitter.addWidget(param_widget)
        splitter.addWidget(self.txt_log)
        splitter.setStretchFactor(1, 3)

        main_layout = QHBoxLayout()
        main_layout.addWidget(splitter)
        self.setLayout(main_layout)

    def _browse_ref(self):
        path, _ = QFileDialog.getOpenFileName(self, "选择 GenBank 文件", "", "GenBank (*.gb *.genbank);;All Files (*)")
        if path:
            self.edit_ref.setText(path)

    def _browse_outdir(self):
        path = QFileDialog.getExistingDirectory(self, "选择输出目录", "")
        if path:
            self.edit_outdir.setText(path)

    def _browse_sample_dir(self):
        path = QFileDialog.getExistingDirectory(self, "选择样本目录", "")
        if path:
            self.edit_sample_dir.setText(path)

    def _on_run(self):
        ref = self.edit_ref.text().strip()
        outdir = self.edit_outdir.text().strip()
        samples = self.edit_sample_dir.text().strip()
        threads = self.spin_threads.value()

        if not (ref and os.path.isfile(ref)):
            QMessageBox.warning(self, "错误", "请选择合法的参考 GenBank 文件")
            return
        if not outdir or not samples:
            QMessageBox.warning(self, "错误", "请完整填写参数")
            return

        os.makedirs(outdir, exist_ok=True)

        script_dir = os.path.dirname(os.path.realpath(__file__))
        core_py = os.path.join(script_dir, "core.py")
        if not os.path.isfile(core_py):
            QMessageBox.critical(self, "错误", f"找不到 core.py：{core_py}")
            return

        cmd = [sys.executable, core_py, "--ref", ref, "--samples", samples, "--outdir", outdir, "--threads", str(threads)]
        self.txt_log.clear()
        self.txt_log.append(">>> 开始执行功能一流程：\n" + " ".join(cmd))

        self.worker = WorkerThread(cmd, workdir=script_dir)
        self.worker.output_signal.connect(self._append_output)
        self.worker.error_signal.connect(self._append_error)
        self.worker.finished_signal.connect(self._proc_finished)
        self.worker.start()

    def _append_output(self, text):
        self.txt_log.moveCursor(QTextCursor.End)
        self.txt_log.insertPlainText(text)

    def _append_error(self, text):
        self.txt_log.moveCursor(QTextCursor.End)
        self.txt_log.setTextColor(QColor("blue"))
        self.txt_log.insertPlainText(text)
        self.txt_log.setTextColor(QColor("black"))

    def _proc_finished(self, exit_code):
        QMessageBox.information(self, "完成", f"功能一执行完成，退出码：{exit_code}")
        self.btn_run.setEnabled(True)

class Function2Tab(QWidget):
    def __init__(self):
        super().__init__()

        self.worker = None
        self._init_ui()

    def _init_ui(self):
        # ===== 左边输入参数区域 =====
        param_layout = QVBoxLayout()

        lbl_complete = QLabel("完整样本目录：")
        self.edit_complete = QLineEdit()
        btn_complete = QPushButton("浏览")
        btn_complete.clicked.connect(self._browse_complete)
        hbox_complete = QHBoxLayout()
        hbox_complete.addWidget(self.edit_complete)
        hbox_complete.addWidget(btn_complete)
        param_layout.addWidget(lbl_complete)
        param_layout.addLayout(hbox_complete)

        lbl_incomplete = QLabel("不完整样本目录：")
        self.edit_incomplete = QLineEdit()
        btn_incomplete = QPushButton("浏览")
        btn_incomplete.clicked.connect(self._browse_incomplete)
        hbox_incomplete = QHBoxLayout()
        hbox_incomplete.addWidget(self.edit_incomplete)
        hbox_incomplete.addWidget(btn_incomplete)
        param_layout.addWidget(lbl_incomplete)
        param_layout.addLayout(hbox_incomplete)

        lbl_output = QLabel("输出目录：")
        self.edit_output = QLineEdit()
        btn_output = QPushButton("浏览")
        btn_output.clicked.connect(self._browse_output)
        hbox_output = QHBoxLayout()
        hbox_output.addWidget(self.edit_output)
        hbox_output.addWidget(btn_output)
        param_layout.addWidget(lbl_output)
        param_layout.addLayout(hbox_output)

        lbl_genes = QLabel("基因列表文件：")
        self.edit_genes = QLineEdit()
        btn_genes = QPushButton("浏览")
        btn_genes.clicked.connect(self._browse_genes)
        hbox_genes = QHBoxLayout()
        hbox_genes.addWidget(self.edit_genes)
        hbox_genes.addWidget(btn_genes)
        param_layout.addWidget(lbl_genes)
        param_layout.addLayout(hbox_genes)

        lbl_threads = QLabel("线程数：")
        self.spin_threads = QSpinBox()
        self.spin_threads.setRange(1, 64)
        self.spin_threads.setValue(4)
        param_layout.addWidget(lbl_threads)
        param_layout.addWidget(self.spin_threads)

        self.chk_skip_extract = QCheckBox("跳过提取")
        self.chk_skip_alignment = QCheckBox("跳过比对")
        self.chk_skip_trimming = QCheckBox("跳过切齐")
        param_layout.addWidget(self.chk_skip_extract)
        param_layout.addWidget(self.chk_skip_alignment)
        param_layout.addWidget(self.chk_skip_trimming)

        lbl_gap = QLabel("Gap阈值：")
        self.spin_gap = QDoubleSpinBox()
        self.spin_gap.setDecimals(2)
        self.spin_gap.setRange(0.0, 1.0)
        self.spin_gap.setValue(0.5)
        param_layout.addWidget(lbl_gap)
        param_layout.addWidget(self.spin_gap)

        self.btn_run = QPushButton("运行")
        self.btn_run.clicked.connect(self._on_run)
        param_layout.addWidget(self.btn_run)
        param_layout.addStretch()

        # ===== 右边日志输出区域 =====
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)

        splitter = QSplitter(Qt.Horizontal)
        param_widget = QWidget()
        param_widget.setLayout(param_layout)
        splitter.addWidget(param_widget)
        splitter.addWidget(self.txt_log)
        splitter.setStretchFactor(1, 3)

        main_layout = QHBoxLayout()
        main_layout.addWidget(splitter)
        self.setLayout(main_layout)

    def _browse_complete(self):
        path = QFileDialog.getExistingDirectory(self, "选择完整样本目录", "")
        if path:
            self.edit_complete.setText(path)

    def _browse_incomplete(self):
        path = QFileDialog.getExistingDirectory(self, "选择不完整样本目录", "")
        if path:
            self.edit_incomplete.setText(path)

    def _browse_output(self):
        path = QFileDialog.getExistingDirectory(self, "选择输出目录", "")
        if path:
            self.edit_output.setText(path)

    def _browse_genes(self):
        path, _ = QFileDialog.getOpenFileName(self, "选择基因列表文件", "", "Text (*.txt);;All Files (*)")
        if path:
            self.edit_genes.setText(path)

    def _on_run(self):
        complete = self.edit_complete.text().strip()
        incomplete = self.edit_incomplete.text().strip()
        output = self.edit_output.text().strip()
        genes = self.edit_genes.text().strip()
        threads = self.spin_threads.value()
        gap = self.spin_gap.value()

        if not (complete and os.path.isdir(complete)):
            QMessageBox.warning(self, "错误", "请选择合法完整样本目录")
            return
        if not (incomplete and os.path.isdir(incomplete)):
            QMessageBox.warning(self, "错误", "请选择合法不完整样本目录")
            return
        if not output or not (genes and os.path.isfile(genes)):
            QMessageBox.warning(self, "错误", "请完整填写参数")
            return

        os.makedirs(output, exist_ok=True)

        script_dir = os.path.dirname(os.path.realpath(__file__))
        cptools_py = os.path.join(script_dir, "CPTOOLS.py")
        if not os.path.isfile(cptools_py):
            QMessageBox.critical(self, "错误", f"找不到 CPTOOLS.py：{cptools_py}")
            return

        cmd = [
            sys.executable, cptools_py,
            "--complete", complete,
            "--incomplete", incomplete,
            "--output", output,
            "--genes", genes,
            "--threads", str(threads),
            "--gap-threshold", str(gap)
        ]

        if self.chk_skip_extract.isChecked():
            cmd.append("--skip-extract")
        if self.chk_skip_alignment.isChecked():
            cmd.append("--skip-alignment")
        if self.chk_skip_trimming.isChecked():
            cmd.append("--skip-trimming")

        self.txt_log.clear()
        self.txt_log.append(">>> 开始执行功能二流程：\n" + " ".join(cmd))

        self.worker = WorkerThread(cmd, workdir=script_dir)
        self.worker.output_signal.connect(self._append_output)
        self.worker.error_signal.connect(self._append_error)
        self.worker.finished_signal.connect(self._proc_finished)
        self.worker.start()

    def _append_output(self, text):
        self.txt_log.moveCursor(QTextCursor.End)
        self.txt_log.insertPlainText(text)

    def _append_error(self, text):
        self.txt_log.moveCursor(QTextCursor.End)
        self.txt_log.setTextColor(QColor("blue"))
        self.txt_log.insertPlainText(text)
        self.txt_log.setTextColor(QColor("black"))

    def _proc_finished(self, exit_code):
        QMessageBox.information(self, "完成", f"功能二执行完成，退出码：{exit_code}")
        self.btn_run.setEnabled(True)

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("叶绿体基因挖掘工具")
        self.resize(1100, 700)

        tabs = QTabWidget()
        tabs.addTab(Function1Tab(), "功能一：测序数据基因挖掘")
        tabs.addTab(Function2Tab(), "功能二：叶绿体组装数据提取")

        layout = QVBoxLayout()
        layout.addWidget(tabs)
        self.setLayout(layout)


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
