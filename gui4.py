#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CPTOOLS GUI — patched version
- 保留原有布局与字段，不擅自删除；仅做对接与修正。
- 模块一对接 core.py（大模块一主控）。
- 模块二对接 CPTOOLS9.py（大模块二主控）。

变更要点：
1) 对“参考 GenBank”和“样本清单 (CSV)”做文件级强校验（必须为单个存在的文件）。
2) 不提供 --species-name 的情况下，必须提供 --ref（否则提前在 GUI 报错）。
3) 修正传参为 --enable_intergenic（与 core.py 一致）。
4) 输出目录可自动创建；日志打印更清晰。
"""

import sys
import os
from PyQt5.QtCore import Qt, QProcess
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QStackedWidget, QSplitter, QGroupBox, QGridLayout,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QTextEdit,
    QFileDialog, QComboBox, QMessageBox
)
from functools import partial

FILE_FILTERS = {
    "参考 GenBank": "GenBank 文件 (*.gb *.gbk *.gbff);;所有文件 (*)",
    "样本清单 (CSV)": "CSV 文件 (*.csv);;所有文件 (*)",
    "构建后报告 (filtered CSV)": "CSV 文件 (*.csv);;所有文件 (*)",
}


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CPMiner v1.0")
        self.resize(1200, 700)

        self.stack = QStackedWidget(self)
        self.stack.addWidget(self._make_welcome_page())
        self.stack.addWidget(self._make_module1_page())
        self.stack.addWidget(self._make_module2_page())

        lay = QVBoxLayout(self)
        lay.addWidget(self.stack)

        self.SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
        self.CORE_SCRIPT = os.path.join(self.SCRIPT_DIR, "core_win.py")
        self.CPT9_SCRIPT = os.path.join(self.SCRIPT_DIR, "CPTOOLS9_win.py")

        self.process = None

    # ---------------- UI PAGES ---------------- #
    def _make_welcome_page(self):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(70, 40, 40, 40)
        v.setSpacing(15)

        title = QLabel("<h1>欢迎使用 CPMiner</h1>")
        title.setAlignment(Qt.AlignCenter)
        v.addWidget(title)

        subtitle = QLabel("版本 1.0")
        subtitle.setAlignment(Qt.AlignCenter)
        v.addWidget(subtitle)

        tip = QLabel("请选择下面的功能模块启动，或点击返回进入其他模块")
        tip.setAlignment(Qt.AlignCenter)
        v.addWidget(tip)

        v.addStretch()

        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(30)
        btn_layout.addStretch()

        btn1 = QPushButton("▶ Module 1：Reads based")
        btn1.setFixedHeight(50)
        btn1.setFixedWidth(240)
        btn1.clicked.connect(lambda: self.stack.setCurrentIndex(1))
        btn_layout.addWidget(btn1)

        btn2 = QPushButton("▶ Module 2：Assembly based")
        btn2.setFixedHeight(50)
        btn2.setFixedWidth(240)
        btn2.clicked.connect(lambda: self.stack.setCurrentIndex(2))
        btn_layout.addWidget(btn2)

        btn_layout.addStretch()
        v.addLayout(btn_layout)

        v.addStretch()
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
            # —— per_gene_filter 参数 ——
            ("最小覆盖深度 (min-depth)", "int", 50),
            ("最大覆盖深度 (max-depth)", "int", 768),
            ("最大数据量Mb (max-size)", "int", 6),
            ("K-mer 大小 (kf)", "int", 31),
            ("组装并行数", "int", 1),
            # —— 其它 ——
            ("迭代次数", "int", 8192),
            ("样本覆盖率", "float", 0.7),
            ("最小比对长度", "int", 500),
            ("最大缺口率%", "float", 60.0),
            ("最小变异位点", "int", 0),
            ("最大变异位点", "int", 200),
            ("最小变异%", "float", 0.0),
            ("最大变异%", "float", 10.0),
            ("最小PI百分比%", "float", 1.0),
            ("构建系统发育树", "check", False),
            ("系统发育树模式", "combo", (["concat", "coalescent"], "concat")),
            ("位点报告(filtered CSV)", "file", None),
            ("建树线程数", "int", 4),  # 仅展示；core.py 的建树线程使用 --threads
        ]
        return self._make_module_page("Module 1：Reads based", fields)

    def _make_module2_page(self):
        fields = [
            ("完整样本目录",       "dir",   None),
            ("不完整样本目录",     "dir",   None),
            ("输出目录",           "dir",   None),
            # BLAST 参数占位
            ("BLAST E-value",      "text",  "1e-5"),
            ("BLAST 身份%",        "float", 70.0),
            ("BLAST 覆盖%",        "float", 70.0),
            # 提取选项
            ("提取间隔区",         "check", False),
            ("软边界 L (bp)",     "int",   200),
            ("软边界 R (bp)",     "int",   200),
            ("模糊匹配模式",       "check", False),
            ("建树线程数",         "int",   4),
            # 过滤参数
            ("样本覆盖率",       "float", 0.70),
            ("最小比对长度", "int",   500),
            ("最大缺口率%",  "float", 60.0),
            ("最小变异位点数", "int", 0),
            ("最大变异位点数", "int", 200),
            ("最小变异率%",   "float", 0.0),
            ("最大变异率%",   "float", 10.0),
            ("最小PI百分比%", "float", 1.0),
            ("构建系统发育树",       "check", False),
            ("系统发育树模式",           "combo", (["concat","coalescent"], "concat")),
            ("位点报告(filtered CSV)",  "file",  None),
        ]
        return self._make_module_page("Module 2：Assembly based", fields)

    def _make_module_page(self, title, fields):
        page = QWidget()
        vlay = QVBoxLayout(page)

        back = QPushButton("← 返回")
        back.setFixedHeight(30)
        back.clicked.connect(lambda: self.stack.setCurrentIndex(0))
        vlay.addWidget(back)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setHandleWidth(5)  # Reduce the gap between parameter and log areas (default is ~5-10px; set to 1 for minimal spacing)

        grp = QGroupBox(f"{title} 参数设置")
        grp.setFixedWidth(700)  # Keep as is
        grid = QGridLayout(grp)
        grid.setVerticalSpacing(5)  # Keep small vertical gap

        setattr(self, f"{title}_widgets", {})
        widgets = getattr(self, f"{title}_widgets")

        for idx, (label, ftype, default) in enumerate(fields):
            row = idx // 2  # Row shared by two fields
            col_pair = (idx % 2) * 2  # 0 for left pair, 2 for right pair
            lbl_col = col_pair
            wgt_col = col_pair + 1

            lbl = QLabel(f"{label}：")
            grid.addWidget(lbl, row, lbl_col)

            if ftype in ("file", "dir"):
                le = QLineEdit()
                le.setFixedWidth(125)  # Keep your preferred width
                btn = QPushButton("…")
                btn.setFixedWidth(32)

                btn.clicked.connect(partial(self._browse, ftype, le, label))
                h = QWidget()
                hh = QHBoxLayout(h)
                hh.setContentsMargins(0, 0, 0, 0)
                hh.addWidget(le)
                hh.addWidget(btn)
                grid.addWidget(h, row, wgt_col)
                widgets[label] = le

            elif ftype == "text":
                le = QLineEdit(str(default) if default is not None else "")
                le.setFixedWidth(125)  # Keep your preferred width
                grid.addWidget(le, row, wgt_col)
                widgets[label] = le

            elif ftype == "int":
                sb = QSpinBox()
                sb.setRange(0, 10**9)
                sb.setValue(int(default or 0))
                sb.setFixedWidth(125)  # Keep your preferred width
                grid.addWidget(sb, row, wgt_col)
                widgets[label] = sb

            elif ftype == "float":
                db = QDoubleSpinBox()
                db.setRange(0, 10**9)
                db.setDecimals(4)
                db.setValue(float(default or 0.0))
                db.setFixedWidth(125)  # Keep your preferred width
                if label.endswith("%"):
                    db.setSuffix(" %")
                grid.addWidget(db, row, wgt_col)
                widgets[label] = db

            elif ftype == "check":
                cb = QCheckBox()
                cb.setChecked(bool(default))
                grid.addWidget(cb, row, wgt_col)
                widgets[label] = cb

            elif ftype == "combo":
                items, sel = default
                cbx = QComboBox()
                cbx.addItems(items)
                cbx.setCurrentText(sel)
                cbx.setFixedWidth(125)  # Keep your preferred width
                grid.addWidget(cbx, row, wgt_col)
                widgets[label] = cbx

        # Add stretch to the row after the last parameter row to push content up
        last_row = (len(fields) - 1) // 2 + 1
        grid.setRowStretch(last_row, 1)

        splitter.addWidget(grp)

        grp_log = QGroupBox("日志输出")
        grp_log.setMinimumWidth(200)  # Keep minimum
        vlog = QVBoxLayout(grp_log)
        te = QTextEdit()
        te.setReadOnly(True)
        vlog.addWidget(te)
        vlog.setStretch(0, 1)  # Ensure QTextEdit fills the group
        splitter.addWidget(grp_log)

        # Increase initial log area width (e.g., 800px; adjust based on your window width of 1000px)
        splitter.setSizes([300, 800])  # [parameter width, log width]; total close to window width minus margins
        splitter.setStretchFactor(0, 0)  # Parameter group does not stretch
        splitter.setStretchFactor(1, 1)  # Log group stretches

        vlay.addWidget(splitter)

        run = QPushButton(f"运行 {title}")
        run.setFixedHeight(40)
        run.clicked.connect(lambda _, t=title, log=te: self.on_run(t, log))
        vlay.addWidget(run)

        return page
    # ---------------- RUNNERS ---------------- #
    def on_run(self, title, log_widget):
        log_widget.clear()
        log_widget.append(f"=== 开始执行 {title} ===")

        widgets = getattr(self, f"{title}_widgets")
        params = {}
        for label, w in widgets.items():
            if isinstance(w, QLineEdit):
                params[label] = w.text()
            elif isinstance(w, (QSpinBox, QDoubleSpinBox)):
                params[label] = w.value()
            elif isinstance(w, QCheckBox):
                params[label] = w.isChecked()
            elif isinstance(w, QComboBox):
                params[label] = w.currentText()
            else:
                params[label] = None
            log_widget.append(f"{label}: {params[label]}")

        if title.startswith("Module 1"):
            # 强校验：必须是“文件”而不是目录
            csv_path = params["样本清单 (CSV)"]
            ref_path = params["参考 GenBank"]
            out_dir  = params["输出目录"]

            if not csv_path or not os.path.isfile(csv_path):
                return self._error("样本清单必须是存在的 CSV 文件（而不是文件夹）")
            if not csv_path.lower().endswith(".csv"):
                return self._error("样本清单文件扩展名应为 .csv")

            if not ref_path or not os.path.isfile(ref_path):
                return self._error("参考 GenBank 必须是存在的文件（而不是文件夹）；当前 GUI 未提供 --species-name，必须选择一个 GenBank 文件")
            if not ref_path.lower().endswith((".gb", ".gbk", ".gbff")):
                return self._error("参考 GenBank 文件扩展名应为 .gb / .gbk / .gbff")

            if not out_dir:
                return self._error("请指定输出目录")
            if not os.path.isdir(out_dir):
                try:
                    os.makedirs(out_dir, exist_ok=True)
                except Exception as e:
                    return self._error(f"无法创建输出目录：{e}")

            self.run_module1(params, log_widget)
        else:
            if not params["完整样本目录"] or not os.path.isdir(params["完整样本目录"]):
                return self._error("请提供存在的“完整样本目录”")
            if not params["不完整样本目录"] or not os.path.isdir(params["不完整样本目录"]):
                return self._error("请提供存在的“不完整样本目录”")
            out_dir = params["输出目录"]
            if not out_dir:
                return self._error("请指定输出目录")
            if not os.path.isdir(out_dir):
                try:
                    os.makedirs(out_dir, exist_ok=True)
                except Exception as e:
                    return self._error(f"无法创建输出目录：{e}")
            self.run_module2(params, log_widget)

    def run_module1(self, p, log):
        args = [
            sys.executable,
            self.CORE_SCRIPT,
            "--file-list", p["样本清单 (CSV)"],
            "--outdir",     p["输出目录"],
            "--threads",    str(p["线程数 (过滤/组装)"]),
        ]
        if p["参考 GenBank"]:
            args += ["--ref", p["参考 GenBank"]]

        if p["提取间隔区"]:
            # 与 core.py 对齐（下划线）
            args.append("--enable_intergenic")
            if p["模糊匹配模式"]:
                args.append("--alias_mode")

        # 软边界
        args += ["-soft_boundary", str(p["软边界 L (bp)"]), str(p["软边界 R (bp)"])]

        # assembler 与 kmer 范围
        args += [
            "-ka",    str(p["kmer (ka)"]),
            "-k_max", str(p["最大 kmer"]),
            "-k_min", str(p["最小 kmer"]),
            "-iteration", str(p["迭代次数"]),
            "--processes-assembler", str(p["组装并行数"]),
        ]

        # per_gene_filter 真实参数
        args += [
            "--min-depth", str(p["最小覆盖深度 (min-depth)"]),
            "--max-depth", str(p["最大覆盖深度 (max-depth)"]),
            "--max-size",  str(p["最大数据量Mb (max-size)"]),
            "-kf",         str(p["K-mer 大小 (kf)"]),
        ]

        # 过滤参数
        args += [
            "--filtered",
            "--seq_pct",       str(p["样本覆盖率"]),
            "--aln_len_min",   str(p["最小比对长度"]),
            "--gap_pct_max",   str(p["最大缺口率%"]),
            "--var_sites_min", str(p["最小变异位点"]),
            "--var_sites_max", str(p["最大变异位点"]),
            "--var_pct_min",   str(p["最小变异%"]),
            "--var_pct_max",   str(p["最大变异%"]),
            "--pi_pct_min",    str(p["最小PI百分比%"]),
        ]

        # 建树
        if p["构建系统发育树"]:
            args += ["--tree", "--tree-mode", p["系统发育树模式"]]
            if p["位点报告(filtered CSV)"]:
                args += ["--report", p["位点报告(filtered CSV)"]]
            # run_tree.py 的线程数使用 --threads；不单独再传

        self._start_process(args, log)

    def run_module2(self, p, log):
        args = [
            sys.executable,
            self.CPT9_SCRIPT,
            "--complete",   p["完整样本目录"],
            "--incomplete", p["不完整样本目录"],
            "--output",     p["输出目录"],
            "--threads",    str(p["建树线程数"]),
        ]
        if p["提取间隔区"]:
            args.append("--enable_intergenic")
        args += ["--soft_boundary", str(p["软边界 L (bp)"]), str(p["软边界 R (bp)"])]
        if p["模糊匹配模式"]:
            args.append("--alias_mode")
        # 过滤参数
        args += [
            "--filtered",
            "--seq_pct",       str(p["样本覆盖率"]),
            "--aln_len_min",   str(p["最小比对长度"]),
            "--gap_pct_max",   str(p["最大缺口率%"]),
            "--var_sites_min", str(p["最小变异位点数"]),
            "--var_sites_max", str(p["最大变异位点数"]),
            "--var_pct_min",   str(p["最小变异率%"]),
            "--var_pct_max",   str(p["最大变异率%"]),
            "--pi_pct_min",    str(p["最小PI百分比%"]),
        ]
        if p["构建系统发育树"]:
            args += ["--tree", "--tree-mode", p["系统发育树模式"]]
            if p["位点报告(filtered CSV)"]:
                args += ["--report", p["位点报告(filtered CSV)"]]

        self._start_process(args, log)
        
    def _browse(self, ft, edit, label_key):
        """统一处理 file/dir 的浏览按钮。"""
        if ft == "file":
            # 根据字段名使用相应的后缀过滤器；找不到就用“所有文件”
            ffilter = FILE_FILTERS.get(label_key, "所有文件 (*)")
            path, _ = QFileDialog.getOpenFileName(self, f"选择 {label_key}", "", ffilter)
        else:
            path = QFileDialog.getExistingDirectory(self, f"选择 {label_key}")
        if path:
            edit.setText(path)

    # ---------------- helpers ---------------- #
    def _start_process(self, args, log_widget):
        cmd_str = " ".join(args)
        log_widget.append("\n执行命令:\n" + cmd_str + "\n")

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyRead.connect(lambda: log_widget.append(bytes(self.process.readAll()).decode('utf-8', errors='ignore')))
        self.process.finished.connect(lambda code, _sig: log_widget.append(f"\n=== 子进程退出：{code} ==="))
        self.process.start(args[0], args[1:])

    def _error(self, msg: str):
        QMessageBox.critical(self, "参数错误", msg)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())
