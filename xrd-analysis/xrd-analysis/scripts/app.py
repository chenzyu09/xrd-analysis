"""Browser interface: drag in XRD files, get the figure, tables and report.

Run:  streamlit run scripts/app.py
"""
import io
import os
import sys
import tempfile
import zipfile

import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xrdtool.pipeline import PALETTE, merge_config, run  # noqa: E402

st.set_page_config(page_title="XRD analysis", layout="wide")
st.title("XRD 分析与论文作图 / XRD analysis")
st.caption("Rigaku .txt/.ras/.rasx/.raw · Bruker .raw · .xy · .xlsx  —  数据只在本机处理 / data stay on this computer")

files = st.file_uploader("XRD 数据文件 / data files", accept_multiple_files=True,
                         type=["txt", "ras", "rasx", "raw", "xy", "xye", "dat", "csv", "xlsx"])
refs = st.file_uploader("参考卡片（CIF 或峰位表，可选）/ references (CIF or stick list, optional)",
                        accept_multiple_files=True, type=["cif", "txt", "csv"])

with st.sidebar:
    st.header("设置 / Settings")
    xr = st.text_input("2θ 范围 / range (e.g. 10 80, empty = auto)", "")
    style = st.selectbox("参考画法 / reference style", ["sticks", "profile", "none"])
    scale = st.selectbox("强度刻度 / intensity scale", ["linear", "sqrt", "log"])
    sub_bg = st.checkbox("扣背景后作图 / plot background-subtracted", False)
    smooth = st.number_input("平滑窗口 (°) / smoothing", 0.0, 1.0, 0.0, 0.01)
    ratio = st.text_input("强度比 / ratio (e.g. 003/104 or 002/101)", "")
    inst = st.number_input("仪器展宽 FWHM (°) / instrument FWHM", 0.0, 0.5, 0.0, 0.01)
    order = st.selectbox("堆叠顺序 / stack order", ["bottom_up", "top_down"])
    font = st.slider("字号 / font size", 8, 24, 14)
    st.header("物相识别 / Phase ID")
    elems = st.text_input("COD 在线检索：样品元素 / elements (e.g. Zn,S,O,H)", "",
                          help="填写后联网到 COD 下载这些元素组成的已知结构并自动比对；留空则用内置库")
    cod_text = st.text_input("COD 关键词（可选）/ keyword (e.g. MIL-101)", "")
    cod_ids = st.text_input("COD 编号（可选，逗号分隔）/ COD IDs", "")

if files:
    labels = []
    cols = st.columns(min(4, len(files)))
    for i, f in enumerate(files):
        with cols[i % len(cols)]:
            labels.append(st.text_input(f"图例 {i + 1}", os.path.splitext(f.name)[0], key=f"lab{i}"))
    if st.button("开始分析 / Run", type="primary"):
        work = tempfile.mkdtemp()
        paths = []
        for f in files:
            p = os.path.join(work, f.name)
            open(p, "wb").write(f.getbuffer())
            paths.append(p)
        cfg = merge_config({})
        cfg["references"] = []
        for r in refs or []:
            p = os.path.join(work, r.name)
            open(p, "wb").write(r.getbuffer())
            cfg["references"].append({"file": p, "label": os.path.splitext(r.name)[0], "show": True})
        if xr.strip():
            try:
                a, b = [float(v) for v in xr.replace(",", " ").split()[:2]]
                cfg["x_range"] = [a, b]
            except ValueError:
                st.warning("范围格式应为: 10 80")
        cfg.update({"reference_style": style, "intensity_scale": scale, "subtract_background_in_plot": sub_bg,
                    "smooth": smooth, "ratios": [ratio] if ratio.strip() else [], "instrument_fwhm": inst,
                    "stack_order": order, "font_size": font})
        sc = cfg["search"]
        if elems.strip():
            sc.update(source="cod", elements=[e for e in elems.replace(" ", ",").replace("，", ",").split(",") if e],
                      text=cod_text.strip(), add_to_plot=True)
        if cod_ids.strip():
            sc["cod_ids"] = [c for c in cod_ids.replace(" ", ",").replace("，", ",").split(",") if c]
        cfg["samples"] = [{"file": p, "label": l, "color": PALETTE[i % len(PALETTE)], "show": True}
                          for i, (p, l) in enumerate(zip(paths, labels))]
        out = os.path.join(work, "xrd_output")
        try:
            with st.spinner("分析中…联网检索时第一次可能需要几分钟 / analysing (first COD search can take minutes)"):
                res = run(paths, cfg, out)
        except Exception as e:  # noqa: BLE001
            st.error(f"分析失败 / failed: {e}")
            st.stop()
        st.image(os.path.join(out, "xrd_stacked.png"))
        for msg in res["results"][0].get("search_log", []) if res["results"] else []:
            if "failed" in msg or "falling back" in msg:
                st.warning(msg)
        for r in res["results"]:
            if r["checks"]:
                st.warning(f"**{r['name']}**: " + " | ".join(r["checks"]))
        st.markdown(open(os.path.join(out, "report.md"), encoding="utf-8").read())
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for root, _, fs in os.walk(out):
                for f in fs:
                    full = os.path.join(root, f)
                    z.write(full, os.path.relpath(full, out))
        st.download_button("下载全部结果 / download all (zip)", buf.getvalue(), "xrd_output.zip")
else:
    st.info("把 XRD 文件拖到上面开始。/ Drop XRD files above to start.")
