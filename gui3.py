#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from PyQt5.QtCore import Qt, QProcess
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QStackedWidget, QSplitter, QGroupBox, QFormLayout, QGridLayout,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QTextEdit,
    QFileDialog, QComboBox
)

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CPTOOLS GUI")
        # 横板 1000×700
        self.resize(1000, 700)

        self.stack = QStackedWidget(self)
        self.stack.addWidget(self._make_welcome_page())
        self.stack.addWidget(self._make_module1_page())
        self.stack.addWidget(self._make_module2_page())

        lay = QVBoxLayout(self)
        lay.addWidget(self.stack)

        self.SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
        self.CORE_SCRIPT = os.path.join(self.SCRIPT_DIR, "core.py")
        self.CPT9_SCRIPT = os.path.join(self.SCRIPT_DIR, "CPTOOLS9.py")

        self.process = None

    def _make_welcome_page(self):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(70, 40, 40, 40)
        v.setSpacing(15)

        # 1) 标题
        title = QLabel("<h1>欢迎使用 CPTOOLS GUI</h1>")
        title.setAlignment(Qt.AlignCenter)
        v.addWidget(title)

        # 2) 子标题 / 版本信息
        subtitle = QLabel("版本 1.0")
        subtitle.setAlignment(Qt.AlignCenter)
        v.addWidget(subtitle)

        # 3) 简短说明
        tip = QLabel("请选择下面的功能模块启动，或点击右下角“帮助”查看使用指南")
        tip.setAlignment(Qt.AlignCenter)
        v.addWidget(tip)

        v.addStretch()

        # 4) 并排按钮
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(30)
        btn_layout.addStretch()

        btn1 = QPushButton("▶ 模块一：测序数据挖掘")
        btn1.setFixedHeight(50)
        btn1.setFixedWidth(240)
        btn1.clicked.connect(lambda: self.stack.setCurrentIndex(1))
        btn_layout.addWidget(btn1)

        btn2 = QPushButton("▶ 模块二：组装数据提取")
        btn2.setFixedHeight(50)
        btn2.setFixedWidth(240)
        btn2.clicked.connect(lambda: self.stack.setCurrentIndex(2))
        btn_layout.addWidget(btn2)

        btn_layout.addStretch()
        v.addLayout(btn_layout)

        v.addStretch()

        # 5) 底部帮助按钮（可选）
        #help_btn = QPushButton("使用帮助")
        #help_btn.setFixedHeight(30)
        #help_btn.clicked.connect(self.show_help)  # 你可以自己实现 show_help()
        #h = QHBoxLayout()
        #h.addStretch()
        #h.addWidget(help_btn)
        #v.addLayout(h)

        return page



    def _make_module1_page(self):
        fields = [
            ("参考 GenBank", "file", None),
            ("样本清单 (CSV)", "file", None),
            ("输出目录", "dir", None),
            ("提取间隔区", "check", False),
            ("软边界 L (bp)", "int", 200),
            ("软边界 R (bp)", "int", 200),
            ("模糊匹配模式", "check", False),
            ("线程数 (过滤/组装)", "int", 4),
            ("kmer (ka)", "int", 39),
            ("最小 kmer", "int", 21),
            ("最大 kmer", "int", 39),
            ("最小深度", "int", 2),
            ("迭代次数", "int", 8192),
            ("样本覆盖率 (seq_pct)", "float", 0.7),
            ("最小比对长度 (aln_len_min)", "int", 500),
            ("最大缺口率% (gap_pct_max)", "float", 60.0),
            ("最小变异位点 (var_sites_min)", "int", 0),
            ("最大变异位点 (var_sites_max)", "int", 200),
            ("最小变异% (var_pct_min)", "float", 0.0),
            ("最大变异% (var_pct_max)", "float", 10.0),
            ("最小 PI 百分比 (pi_pct_min)", "float", 1.0),
            ("运行后构建系统发育树", "check", False),
            ("系统发育树模式", "combo", (["concat", "coalescent"], "concat")),
            ("构建后报告 (filtered CSV)", "file", None),
            ("树构建线程数", "int", 4),
        ]
        return self._make_module_page("模块一：测序数据挖掘", fields)

    def _make_module2_page(self):
        fields=[
            # I/O 路径
            ("完整样本目录",       "dir",   None),
            ("不完整样本目录",     "dir",   None),
            ("输出目录",           "dir",   None),
            # BLAST 参数
            ("BLAST E-value",      "text",  "1e-5"),
            ("BLAST 身份%",        "float", 70.0),
            ("BLAST 覆盖%",        "float", 70.0),
            # 提取选项
            ("提取间隔区",         "check", False),
            ("软边界 L (bp)",     "int",   200),
            ("软边界 R (bp)",     "int",   200),
            ("模糊匹配模式",       "check", False),
            ("建树线程数",         "int",   4),
            # —— 新增：过滤 summary_table 的参数 —— #
            ("样本覆盖率 (seq_pct)",      "float", 0.70),
            ("最小比对长度 (aln_len_min)","int",   500),
            ("最大缺口率% (gap_pct_max)", "float", 60.0),
            ("最小变异位点数 (var_sites_min)","int",   0),
            ("最大变异位点数 (var_sites_max)","int",   200),
            ("最小变异率% (var_pct_min)",  "float", 0.0),
            ("最大变异率% (var_pct_max)",  "float", 10.0),
            ("最小 PI 百分比 (pi_pct_min)","float", 1.0),
            # —— 新增：建树参数 —— #
            ("运行后构建系统发育树",    "check", False),
            ("系统发育树模式",        "combo", (["concat","coalescent"], "concat")),
            ("构建后报告 (filtered CSV)","file",  None),
            ("树构建线程数",          "int",   4),
        ]
        return self._make_module_page("模块二：组装数据提取", fields)



    def _make_module_page(self, title, fields):
        page = QWidget()
        vlay = QVBoxLayout(page)

        # —— 返回按钮 —— #
        back = QPushButton("← 返回")
        back.setFixedHeight(30)
        back.clicked.connect(lambda: self.stack.setCurrentIndex(0))
        vlay.addWidget(back)

        splitter = QSplitter(Qt.Vertical)

        # —— 参数区 —— #
        grp = QGroupBox(f"{title} 参数设置")
        grid = QGridLayout(grp)
        # 四列： [标签1, 控件1, 标签2, 控件2]
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)

        setattr(self, f"{title}_widgets", {})
        widgets = getattr(self, f"{title}_widgets")

        n_rows = (len(fields) + 1) // 2
        for idx, (label, ftype, default) in enumerate(fields):
            row = idx // 2
            col_pair = idx % 2  # 0=左列, 1=右列
            lbl_col = col_pair * 2
            wgt_col = lbl_col + 1

            # 1) 标签
            lbl = QLabel(f"{label}：")
            grid.addWidget(lbl, row, lbl_col)

            # 2) 控件
            if ftype in ("file", "dir"):
                le = QLineEdit(); le.setMaximumWidth(200)
                btn = QPushButton("…"); btn.setFixedWidth(30)
                def _browse(ft=ftype, edit=le, text=label):
                    if ft=="file":
                        p,_=QFileDialog.getOpenFileName(self, f"选择 {text}")
                    else:
                        p=QFileDialog.getExistingDirectory(self, f"选择 {text}")
                    if p: edit.setText(p)
                btn.clicked.connect(_browse)
                h = QWidget()
                hh = QHBoxLayout(h); hh.setContentsMargins(0,0,0,0)
                hh.addWidget(le); hh.addWidget(btn)
                grid.addWidget(h, row, wgt_col)
                widgets[label] = le

            elif ftype=="text":
                le = QLineEdit(str(default)); le.setMaximumWidth(200)
                grid.addWidget(le, row, wgt_col)
                widgets[label] = le

            elif ftype=="int":
                sb = QSpinBox(); sb.setRange(0,10**9); sb.setValue(default or 0)
                sb.setMaximumWidth(100)
                grid.addWidget(sb, row, wgt_col)
                widgets[label] = sb

            elif ftype=="float":
                db = QDoubleSpinBox(); db.setRange(0,10**9); db.setDecimals(4)
                db.setValue(default or 0.0); db.setMaximumWidth(120)
                if label.endswith("%"):
                    db.setSuffix(" %")
                grid.addWidget(db, row, wgt_col)
                widgets[label] = db

            elif ftype=="check":
                cb = QCheckBox(); cb.setChecked(bool(default))
                grid.addWidget(cb, row, wgt_col)
                widgets[label] = cb

            elif ftype=="combo":
                items, sel = default
                cbx = QComboBox(); cbx.addItems(items)
                cbx.setCurrentText(sel)
                cbx.setMaximumWidth(120)
                grid.addWidget(cbx, row, wgt_col)
                widgets[label] = cbx

        splitter.addWidget(grp)
        # —— 调整 参数区 vs 日志区 的占比 —— #
        # 参数区(第0个子控件) 占 3 份，日志区(第1个子控件) 占 1 份
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        # —— 如果你想要用固定像素高度，也可以： splitter.setSizes([500, 150])
        # —— 日志区 —— #
        grp_log = QGroupBox("日志输出")
        vlog = QVBoxLayout(grp_log)
        te = QTextEdit(); te.setReadOnly(True)
        vlog.addWidget(te)
        splitter.addWidget(grp_log)

        vlay.addWidget(splitter)

        # —— 底部运行按钮 —— #
        run = QPushButton(f"运行 {title}")
        run.setFixedHeight(40)
        run.clicked.connect(lambda _, t=title, log=te: self.on_run(t, log))
        vlay.addWidget(run)

        return page


    def on_run(self, title, log_widget):
        # 1) 先把日志清空
        log_widget.clear()
        log_widget.append(f"=== 开始执行 {title} ===")

        # 2) 收集所有参数到一个 dict
        widgets = getattr(self, f"{title}_widgets")
        params = {}
        for label, w in widgets.items():
            if isinstance(w, QLineEdit):
                params[label] = w.text()
            elif isinstance(w, QSpinBox) or isinstance(w, QDoubleSpinBox):
                params[label] = w.value()
            elif isinstance(w, QCheckBox):
                params[label] = w.isChecked()
            elif isinstance(w, QComboBox):
                params[label] = w.currentText()
            else:
                params[label] = None
            log_widget.append(f"{label}: {params[label]}")

        log_widget.append("… 正在启动子进程 …")

        # 3) 根据标题调用不同脚本
        if title.startswith("模块一"):
            self.run_module1(params, log_widget)
        else:
            self.run_module2(params, log_widget)  

    def run_module1(self, p, log):
        # 构造命令行
        args = [
            sys.executable,
            self.CORE_SCRIPT,
            "--file-list", p["样本清单 (CSV)"],
            "--outdir",     p["输出目录"],
        ]
        # 可选 ref
        if p["参考 GenBank"]:
            args += ["--ref", p["参考 GenBank"]]
        # 提取 intergenic
        if p["提取间隔区"]:
            args.append("--enable-intergenic")
            if p["模糊匹配模式"]:
                args.append("--alias_mode")
        # 软边界
        args += ["-soft_boundary", str(p["软边界 L (bp)"]), str(p["软边界 R (bp)"])]
        # 核心线程
        args += ["--threads", str(p["线程数 (过滤/组装)"])]
        # assembler 参数
        args += [
            "-ka",      str(p["kmer (ka)"]),
            "-k_max",   str(p["最大 kmer"]),
            "-k_min",   str(p["最小 kmer"]),
            "-limit_count", str(p["最小深度"]),
            "-iteration",   str(p["迭代次数"]),
            "--trim-mode",  "automated1",
        ]
        # 过滤 summary_table 参数
        args += [
            "--filtered",
            "--seq_pct",     str(p["样本覆盖率 (seq_pct)"]),
            "--aln_len_min", str(p["最小比对长度 (aln_len_min)"]),
            "--gap_pct_max", str(p["最大缺口率% (gap_pct_max)"]),
            "--var_sites_min", str(p["最小变异位点 (var_sites_min)"]),
            "--var_sites_max", str(p["最大变异位点 (var_sites_max)"]),
            "--var_pct_min",  str(p["最小变异% (var_pct_min)"]),
            "--var_pct_max",  str(p["最大变异% (var_pct_max)"]),
            "--pi_pct_min",   str(p["最小 PI 百分比 (pi_pct_min)"]),
        ]
        # 建树
        if p["运行后构建系统发育树"]:
            args += ["--tree", "--tree-mode", p["系统发育树模式"]]
            if p["构建后报告 (filtered CSV)"]:
                args += ["--report", p["构建后报告 (filtered CSV)"]]
        args += ["-t", str(p["树构建线程数"])]

        log.append("执行命令:\n" + " ".join(args))

        # 启动 QProcess
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyRead.connect(lambda:
            log.append(bytes(self.process.readAll()).decode('utf-8'))
        )
        self.process.start(args[0], args[1:])

    # —— 新增方法：运行功能二 CPTOOLS9.py —— #
    def run_module2(self, p, log):
        args = [
            sys.executable,
            self.CPT9_SCRIPT,
            "--complete",   p["完整样本目录"],
            "--incomplete", p["不完整样本目录"],
            "--output",     p["输出目录"],
        ]
        # BLAST 就用默认 text / float 控件
        # 提取 intergenic
        if p["提取间隔区"]:
            args.append("--enable-intergenic")
        # 软边界
        args += ["--soft_boundary", str(p["软边界 L (bp)"]), str(p["软边界 R (bp)"])]
        if p["模糊匹配模式"]:
            args.append("--alias_mode")
        # 过滤参数同模块一
        args += [
            "--filtered",
            "--seq_pct",     str(p["样本覆盖率 (seq_pct)"]),
            "--aln_len_min", str(p["最小比对长度 (aln_len_min)"]),
            "--gap_pct_max", str(p["最大缺口率% (gap_pct_max)"]),
            "--var_sites_min", str(p["最小变异位点数 (var_sites_min)"]),
            "--var_sites_max", str(p["最大变异位点数 (var_sites_max)"]),
            "--var_pct_min",  str(p["最小变异率% (var_pct_min)"]),
            "--var_pct_max",  str(p["最大变异率% (var_pct_max)"]),
            "--pi_pct_min",   str(p["最小 PI 百分比 (pi_pct_min)"]),
        ]
        # 建树
        if p["运行后构建系统发育树"]:
            args += ["--tree", "--tree-mode", p["系统发育树模式"]]
            if p["构建后报告 (filtered CSV)"]:
                args += ["--report", p["构建后报告 (filtered CSV)"]]
        args += ["--threads", str(p["树构建线程数"])]

        log.append("执行命令:\n" + " ".join(args))

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyRead.connect(lambda:
            log.append(bytes(self.process.readAll()).decode('utf-8'))
        )
        self.process.start(args[0], args[1:])


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())
