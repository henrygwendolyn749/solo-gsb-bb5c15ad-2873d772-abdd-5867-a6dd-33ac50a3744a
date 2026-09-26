"""app.digest 的 /erase-sites 纯标准库测试：python3 tests/test_erase_sites.py"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.digest import DigestError, run_erase_sites  # noqa: E402

PASS = 0
FAIL = 0


def check(name, cond, extra=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name} {extra}")


def expect_error(name, payload, snippet):
    try:
        run_erase_sites(payload)
    except DigestError as e:
        check(name, snippet in e.reason, f"reason={e.reason!r}")
        return
    check(name, False, "未抛出 DigestError")


AAA = {"name": "E", "recognition": "AAA", "cut_offset": 0}
ECORI = {"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}


def payload(**kw):
    base = {
        "sequence": "AAGAAA",
        "topology": "linear",
        "enzymes": [AAA],
        "edits": [{"position": 4, "base": "C"}],
    }
    base.update(kw)
    return base


# --- 两个 EcoRI 位点，两个候选替换缺一不可；回文酶同一起点正反链只计一个位点 ---
r = run_erase_sites({
    "sequence": "AAGAATTCTTTGAATTCGG",
    "topology": "linear",
    "enzymes": [ECORI],
    "edits": [{"position": 3, "base": "T"}, {"position": 12, "base": "T"}],
})
check("ecoRI two-site status", r["status"] == "solution", str(r))
check("ecoRI two-site indices", r["edit_indices"] == [0, 1], str(r))
check("ecoRI two-site edit_count", r["edit_count"] == 2)
check("ecoRI two-site edited", r["edited_sequence"] == "AAGTATTCTTTGTATTCGG",
      r["edited_sequence"])
check("ecoRI enzyme entry", r["enzymes"] == [
    {"name": "EcoRI", "recognition": "GAATTC",
     "original_sites": 2, "remaining_sites": 0}], str(r["enzymes"]))

# --- 原本就没有位点：空替换组即最优（改动数 0） ---
r = run_erase_sites(payload(sequence="ACGTACGT", edits=[]))
check("already clean solution", r["status"] == "solution", str(r))
check("already clean empty indices", r["edit_indices"] == [] and r["edit_count"] == 0,
      str(r))
check("already clean sequence unchanged", r["edited_sequence"] == "ACGTACGT")
check("already clean counts", r["enzymes"][0]["original_sites"] == 0
      and r["enzymes"][0]["remaining_sites"] == 0, str(r))

# --- 单个替换即可：不选更大的组合；多个单替换可行时取下标最小者 ---
r = run_erase_sites(payload(
    sequence="AAA",
    edits=[{"position": 0, "base": "C"}, {"position": 1, "base": "C"},
           {"position": 2, "base": "C"}],
))
check("single edit suffices", r["status"] == "solution"
      and r["edit_indices"] == [0] and r["edit_count"] == 1, str(r))
check("single edit edited seq", r["edited_sequence"] == "CAA", str(r))

# --- 同时应用后重新扫描：消除旧位点却生成新位点的组合被否决 ---
# AAGAAA：位点在 3；改位置 2 (G->A) 得 AAAAAA 反而遍地位点；改位置 4 才干净
r = run_erase_sites(payload(
    edits=[{"position": 2, "base": "A"}, {"position": 4, "base": "C"}],
))
check("rescan rejects site-creating edit", r["status"] == "solution"
      and r["edit_indices"] == [1], str(r))
check("rescan edited seq", r["edited_sequence"] == "AAGACA", str(r))
check("rescan original count", r["enzymes"][0]["original_sites"] == 1, str(r))

# --- 需要两个替换；多个可行二元组中取输入下标字典序最小者 ---
# AAATAAA：位点在 0 与 4；idx0 只消位点4，idx1/idx2 只消位点0，
# 可行对为 {0,1} 与 {0,2}，应选 [0,1]
r = run_erase_sites(payload(
    sequence="AAATAAA",
    edits=[{"position": 5, "base": "C"}, {"position": 1, "base": "C"},
           {"position": 2, "base": "C"}],
))
check("pair needed, lexicographic", r["status"] == "solution"
      and r["edit_indices"] == [0, 1] and r["edit_count"] == 2, str(r))
check("pair edited seq", r["edited_sequence"] == "ACATACA", r["edited_sequence"])
check("pair original count", r["enzymes"][0]["original_sites"] == 2, str(r))

# --- 回文识别序列：同一起点正反链命中计为一个位点 ---
r = run_erase_sites(payload(
    sequence="GAATTC", enzymes=[ECORI],
    edits=[{"position": 0, "base": "A"}],
))
check("palindrome counts once", r["enzymes"][0]["original_sites"] == 1, str(r))
check("palindrome erased", r["status"] == "solution"
      and r["edited_sequence"] == "AAATTC"
      and r["enzymes"][0]["remaining_sites"] == 0, str(r))

# --- 反向互补链上的位点也要计入并消除（AAA 的反向互补为 TTT） ---
r = run_erase_sites(payload(
    sequence="AAATTT",
    edits=[{"position": 1, "base": "G"}, {"position": 4, "base": "A"}],
))
check("reverse-strand site counted", r["enzymes"][0]["original_sites"] == 2, str(r))
check("reverse-strand erased", r["status"] == "solution"
      and r["edit_indices"] == [0, 1]
      and r["edited_sequence"] == "AGATAT", str(r))

# --- 环状：跨首尾（原点）的位点可识别并消除 ---
# 环状 AATAA：AAA 跨原点命中起点 3 与 4，共 2 个位点；改位置 0 后全清
r = run_erase_sites(payload(
    sequence="AATAA", topology="circular",
    edits=[{"position": 0, "base": "C"}],
))
check("circular wrap sites counted", r["enzymes"][0]["original_sites"] == 2, str(r))
check("circular wrap erased", r["status"] == "solution"
      and r["edit_indices"] == [0] and r["edited_sequence"] == "CATAA"
      and r["enzymes"][0]["remaining_sites"] == 0, str(r))

# --- 多酶：所有酶的位点都必须清零，逐酶返回原有/剩余数 ---
r = run_erase_sites(payload(
    sequence="AAACCC",
    enzymes=[{"name": "E1", "recognition": "AAA", "cut_offset": 0},
             {"name": "E2", "recognition": "CCC", "cut_offset": 0}],
    edits=[{"position": 1, "base": "T"}, {"position": 4, "base": "T"}],
))
check("multi-enzyme status", r["status"] == "solution"
      and r["edit_indices"] == [0, 1], str(r))
check("multi-enzyme edited", r["edited_sequence"] == "ATACTC", r["edited_sequence"])
check("multi-enzyme counts", [
    (e["name"], e["original_sites"], e["remaining_sites"]) for e in r["enzymes"]
] == [("E1", 1, 0), ("E2", 1, 0)], str(r["enzymes"]))

# --- 无法清除：候选替换覆盖不到全部位点 ---
r = run_erase_sites(payload(sequence="AAAAAA",
                            edits=[{"position": 0, "base": "C"}]))
check("unsatisfiable no_solution", r == {"status": "no_solution"}, str(r))
r = run_erase_sites(payload(edits=[]))
check("empty edits with sites no_solution", r == {"status": "no_solution"}, str(r))

# --- 大小写归一化：序列与替换碱基均可小写 ---
r = run_erase_sites(payload(sequence="aagaaa",
                            edits=[{"position": 4, "base": "c"}]))
check("lowercase normalized", r["status"] == "solution"
      and r["edited_sequence"] == "AAGACA", str(r))

# --- 边界：500 碱基接受、501 拒绝；8 个酶接受、9 个拒绝；10 个替换接受 ---
r = run_erase_sites(payload(sequence="A" * 500,
                            edits=[{"position": 0, "base": "C"}]))
check("500-length accepted", r == {"status": "no_solution"}, str(r)[:120])
expect_error("501-length rejected", payload(sequence="A" * 501), "上限 500")

enz8 = [dict(AAA, name=f"E{i}") for i in range(8)]
r = run_erase_sites(payload(enzymes=enz8))
check("8 enzymes accepted", r["status"] == "solution"
      and len(r["enzymes"]) == 8, str(r)[:200])
expect_error("9 enzymes rejected", payload(enzymes=enz8 + [AAA]), "酶数量")
expect_error("0 enzymes rejected", payload(enzymes=[]), "酶数量")

edits10 = [{"position": i, "base": "C"} for i in range(10)]
r = run_erase_sites(payload(sequence="A" * 10, edits=edits10))
check("10 edits accepted", r["status"] in ("solution", "no_solution"), str(r)[:120])
expect_error("11 edits rejected",
             payload(sequence="A" * 11,
                     edits=edits10 + [{"position": 10, "base": "C"}]),
             "不能超过 10")

# --- 非法输入：全部抛出 DigestError（HTTP 层映射为 400） ---
expect_error("payload not dict", ["not", "dict"], "JSON 对象")
expect_error("missing sequence", payload(sequence=None), "sequence")
expect_error("illegal base in sequence", payload(sequence="AANAAA"), "非法碱基")
expect_error("empty sequence", payload(sequence="  "), "不能为空")
expect_error("bad topology", payload(topology="plasmid"), "拓扑")
expect_error("missing topology", payload(topology=None), "topology")
expect_error("missing enzymes", payload(enzymes=None), "enzymes")
expect_error("enzyme bad recognition",
             payload(enzymes=[{"name": "E", "recognition": "AAN",
                               "cut_offset": 0}]), "非法碱基")
expect_error("enzyme bad offset",
             payload(enzymes=[{"name": "E", "recognition": "AAA",
                               "cut_offset": 9}]), "越界")
expect_error("missing edits", payload(edits=None), "edits")
expect_error("edits not list", payload(edits="x"), "数组")
expect_error("edit not object", payload(edits=[1]), "对象")
expect_error("position missing", payload(edits=[{"base": "C"}]), "position")
expect_error("position bool",
             payload(edits=[{"position": True, "base": "C"}]), "整数")
expect_error("position float",
             payload(edits=[{"position": 1.5, "base": "C"}]), "整数")
expect_error("position negative",
             payload(edits=[{"position": -1, "base": "C"}]), "越界")
expect_error("position == length",
             payload(edits=[{"position": 6, "base": "C"}]), "越界")
expect_error("duplicate position",
             payload(edits=[{"position": 4, "base": "C"},
                            {"position": 4, "base": "T"}]), "同一位置")
expect_error("base missing", payload(edits=[{"position": 4}]), "base")
expect_error("base not string",
             payload(edits=[{"position": 4, "base": 5}]), "字符串")
expect_error("base illegal",
             payload(edits=[{"position": 4, "base": "N"}]), "A/C/G/T")
expect_error("base too long",
             payload(edits=[{"position": 4, "base": "AC"}]), "单个")
expect_error("base same as original",
             payload(edits=[{"position": 4, "base": "a"}]), "改变原碱基")

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
