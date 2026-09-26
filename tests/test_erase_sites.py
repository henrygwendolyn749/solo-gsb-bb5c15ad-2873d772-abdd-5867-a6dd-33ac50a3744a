"""app.digest 的 /erase-sites 纯标准库测试：python3 tests/test_erase_sites.py"""

import sys
import os
import random
import itertools

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.digest import (  # noqa: E402
    DigestError,
    reverse_complement,
    run_erase_sites,
)

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


ECORI = {"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}
BSAI = {"name": "BsaI", "recognition": "GGTCTC", "cut_offset": 1}


def reference_sites(seq, recognitions, topology):
    """独立的朴素窗口扫描，作为 _find_site_starts 的交叉校验参考。"""
    n = len(seq)
    table = []
    for rec in recognitions:
        L = len(rec)
        rc = reverse_complement(rec)
        search = seq + seq[: L - 1] if topology == "circular" else seq
        starts = range(n) if topology == "circular" else range(max(0, n - L + 1))
        sites = set()
        for s in starts:
            window = search[s : s + L]
            if window == rec or window == rc:
                sites.add(s)
        table.append(sites)
    return table


def apply_edits(seq, combo, edits):
    chars = list(seq)
    for idx in combo:
        chars[edits[idx]["position"]] = edits[idx]["base"]
    return "".join(chars)


def brute_force_optimal(seq, recognitions, topology, edits):
    """与实现相同的择优规则（最少改动、同数下标字典序最小）的朴素枚举参考。"""
    m = len(edits)
    for size in range(0, m + 1):
        for combo in itertools.combinations(range(m), size):
            candidate = apply_edits(seq, combo, edits)
            if not any(reference_sites(candidate, recognitions, topology)):
                return list(combo), candidate
    return None, None


# --- 基本：单酶单位点，单位基替换清除 ---
# TTTGAATTCTTT：GAATTC 起点 3；position 3 为 G，改成 C
r = run_erase_sites({
    "sequence": "TTTGAATTCTTT",
    "topology": "linear",
    "enzymes": [ECORI],
    "edits": [{"position": 3, "base": "C"}],
})
check("basic status", r["status"] == "solution", str(r))
check("basic sequence", r["sequence"] == "TTTCAATTCTTT", str(r["sequence"]))
check("basic indices", r["edit_indices"] == [0] and r["edit_count"] == 1,
      str(r["edit_indices"]))
check("basic original/remaining",
      r["enzyme_sites"][0]["original_sites"] == 1
      and r["enzyme_sites"][0]["remaining_sites"] == 0
      and r["enzyme_sites"][0]["recognition"] == "GAATTC"
      and r["enzyme_sites"][0]["name"] == "EcoRI",
      str(r["enzyme_sites"]))
check("basic length preserved", r["length"] == 12 and len(r["sequence"]) == 12)

# --- 反链（反向互补）命中也算位点，替换可清除 ---
# AAGAGACCAAA：GAGACC（GGTCTC 的反向互补）起点 2
r = run_erase_sites({
    "sequence": "AAGAGACCAAA",
    "topology": "linear",
    "enzymes": [BSAI],
    "edits": [{"position": 2, "base": "T"}],
})
check("reverse strand hit erased", r["status"] == "solution"
      and r["enzyme_sites"][0]["original_sites"] == 1
      and r["enzyme_sites"][0]["remaining_sites"] == 0, str(r))

# --- 回文酶同一起点正反两条链命中，计为一个位点 ---
r = run_erase_sites({
    "sequence": "AAGAATTCAA",
    "topology": "linear",
    "enzymes": [ECORI],
    "edits": [{"position": 2, "base": "T"}],
})
check("palindrome counts as one site",
      r["enzyme_sites"][0]["original_sites"] == 1
      and r["enzyme_sites"][0]["remaining_sites"] == 0, str(r["enzyme_sites"]))

# --- 环状跨原点位点 ---
# TCTCAAAAGG（len 10）：末尾 GG + 开头 TCTC 构成跨原点的 GGTCTC，start=8
r = run_erase_sites({
    "sequence": "TCTCAAAAGG",
    "topology": "circular",
    "enzymes": [BSAI],
    "edits": [
        {"position": 8, "base": "A"},  # 破坏跨原点位点
        {"position": 0, "base": "G"},  # 同样破坏，但下标更大
    ],
})
check("circular wrap site recognized and erased",
      r["status"] == "solution" and r["edit_indices"] == [0]
      and r["enzyme_sites"][0]["original_sites"] == 1
      and r["enzyme_sites"][0]["remaining_sites"] == 0, str(r))
# 同序列线性时没有跨原点位点 -> 空编辑即解
rl = run_erase_sites({
    "sequence": "TCTCAAAAGG",
    "topology": "linear",
    "enzymes": [BSAI],
    "edits": [{"position": 8, "base": "A"}],
})
check("linear no site -> empty edit solution", rl["status"] == "solution"
      and rl["edit_indices"] == [] and rl["edit_count"] == 0
      and rl["sequence"] == "TCTCAAAAGG"
      and rl["enzyme_sites"][0]["original_sites"] == 0, str(rl))

# --- 编辑同时应用后重扫：局部消除却在邻位生成新位点时，该编辑不可行 ---
# BsaI：GGTCTC / rc GAGACC。序列 GGTCTCAGACC（len 11），仅有正向位点 GGTCTC@0；
# 位置 6..10 预置为 AGACC。诱饵编辑 position 5 C->G：@0 失配（GGTCTG），
# 却在 @5 补出反向互补位点 GAGACC（seq[5..10] = GAGACC）。
bait_seq = "GGTCTCAGACC"
check("bait seq length", len(bait_seq) == 11)
pre = run_erase_sites({
    "sequence": bait_seq, "topology": "linear", "enzymes": [BSAI],
    "edits": [{"position": 2, "base": "A"}],  # 干净编辑，借其返回查原始计数
})
check("bait setup: one site @0",
      pre["status"] == "solution"
      and pre["enzyme_sites"][0]["original_sites"] == 1, str(pre))
bait = run_erase_sites({
    "sequence": bait_seq, "topology": "linear", "enzymes": [BSAI],
    "edits": [{"position": 5, "base": "G"}],
})
check("bait alone fails (new reverse-strand site at 5)",
      bait == {"status": "no_solution"}, str(bait))
# 干净编辑 position 2 T->A：GGTCTCAGACC -> GGACTCAGACC，无任一模式
clean = run_erase_sites({
    "sequence": bait_seq, "topology": "linear", "enzymes": [BSAI],
    "edits": [{"position": 2, "base": "A"}],
})
check("clean edit alone solves", clean["status"] == "solution"
      and clean["sequence"] == "GGACTCAGACC"
      and clean["enzyme_sites"][0]["remaining_sites"] == 0, str(clean))
# 两个候选都给出时，诱饵（下标 0）不可行，选下标 1
both = run_erase_sites({
    "sequence": bait_seq, "topology": "linear", "enzymes": [BSAI],
    "edits": [{"position": 5, "base": "G"}, {"position": 2, "base": "A"}],
})
check("rescan rejects bait and picks clean edit",
      both["status"] == "solution" and both["edit_indices"] == [1], str(both))

# --- 多个位点需要多个编辑；单编辑不足时 no_solution ---
# GAGACCGGTCTC：反向位点 @0 与正向位点 @7
multi = "GAGACCGGTCTC"
r = run_erase_sites({
    "sequence": multi, "topology": "linear", "enzymes": [BSAI],
    "edits": [
        {"position": 0, "base": "C"},   # 只消反向链位点 @0
        {"position": 7, "base": "A"},   # 只消正向位点 @7
    ],
})
check("two sites need two edits", r["status"] == "solution"
      and r["edit_indices"] == [0, 1] and r["edit_count"] == 2
      and r["enzyme_sites"][0]["original_sites"] == 2
      and r["enzyme_sites"][0]["remaining_sites"] == 0, str(r))
r_one = run_erase_sites({
    "sequence": multi, "topology": "linear", "enzymes": [BSAI],
    "edits": [{"position": 0, "base": "C"}],
})
check("one edit insufficient -> no_solution",
      r_one == {"status": "no_solution"}, str(r_one))

# --- 同数时下标列表字典序最小 ---
# TTTGAATTCTTT 位点 @3；position 3 与 position 4 的替换各自都能破坏位点
for order in ([3, 4], [4, 3]):
    r = run_erase_sites({
        "sequence": "TTTGAATTCTTT", "topology": "linear", "enzymes": [ECORI],
        "edits": [{"position": p, "base": "C"} for p in order],
    })
    expected_idx = [0]  # 无论位置如何，都选输入下标 0
    check(f"lexicographic singleton order {order}",
          r["status"] == "solution" and r["edit_indices"] == expected_idx,
          str(r))

# 二元组合字典序：位点 @0、@9（GAATTCAAA...GAATTC），编辑 0/1 各自只破坏 @0，
# 编辑 2 只破坏 @9 -> 可行二元组 {0,2} 与 {1,2} 中选 [0,2]
two_site = "GAATTC" + "AAA" + "GAATTC"
r = run_erase_sites({
    "sequence": two_site, "topology": "linear", "enzymes": [ECORI],
    "edits": [
        {"position": 0, "base": "T"},
        {"position": 1, "base": "T"},
        {"position": 9, "base": "T"},
    ],
})
check("two-site lexicographic pair [0,2]", r["status"] == "solution"
      and r["edit_indices"] == [0, 2]
      and r["enzyme_sites"][0]["original_sites"] == 2
      and r["enzyme_sites"][0]["remaining_sites"] == 0, str(r))

# --- 改动数最少优先于下标字典序 ---
# 序列 CGAATTC（len 7）：酶1 GAATTC@1，酶2 CGAATT@0（rc=AATTCG，均不命中反链）。
# 编辑 0（C->T）只消酶2，编辑 1（末位 C->A）只消酶1，二者必须合用（字典序最小
# 的二元解 [0,1]）；编辑 2（G@1->C）使两位点同时失配且不产生新位点 -> 单编辑 [2]
# 规模更小，应优先选中。
few_enzymes = [
    {"name": "E1", "recognition": "GAATTC", "cut_offset": 1},
    {"name": "E2", "recognition": "CGAATT", "cut_offset": 0},
]
check("fewest setup sites",
      reference_sites("CGAATTC", ["GAATTC", "CGAATT"], "linear") == [{1}, {0}])
r = run_erase_sites({
    "sequence": "CGAATTC", "topology": "linear", "enzymes": few_enzymes,
    "edits": [
        {"position": 0, "base": "T"},  # 只消酶2
        {"position": 6, "base": "A"},  # 只消酶1（GAATTC@1 末位变 A）
        {"position": 1, "base": "C"},  # 单编辑同时消两者
    ],
})
check("fewest edits preferred over lexicographic pair",
      r["status"] == "solution" and r["edit_indices"] == [2]
      and r["sequence"] == "CCAATTC"
      and [e["original_sites"] for e in r["enzyme_sites"]] == [1, 1]
      and [e["remaining_sites"] for e in r["enzyme_sites"]] == [0, 0], str(r))
# 对照：去掉单编辑竞争者后，[0,1] 是最小可行组合
r = run_erase_sites({
    "sequence": "CGAATTC", "topology": "linear", "enzymes": few_enzymes,
    "edits": [{"position": 0, "base": "T"}, {"position": 6, "base": "A"}],
})
check("pair [0,1] feasible without singleton competitor",
      r["status"] == "solution" and r["edit_indices"] == [0, 1]
      and r["sequence"] == "TGAATTA", str(r))

# --- 同时应用：多编辑结果等于一次性替换 ---
r = run_erase_sites({
    "sequence": two_site, "topology": "linear", "enzymes": [ECORI],
    "edits": [
        {"position": 0, "base": "T"},
        {"position": 9, "base": "T"},
    ],
})
expected = "T" + two_site[1:9] + "T" + two_site[10:]
check("simultaneous application", r["status"] == "solution"
      and r["sequence"] == expected and r["edit_indices"] == [0, 1],
      str((r["sequence"], expected)))

# --- 多酶：逐酶报告 original/remaining，顺序与输入一致 ---
# GAATTCGGG（len 9）：EcoRI@0、GGG@6；CAT（rc=ATG）始终不命中。
# position 3 A->C 破坏 EcoRI；position 8 G->T 把 GGG@6 变 GGT。
r = run_erase_sites({
    "sequence": "GAATTCGGG",
    "topology": "linear",
    "enzymes": [
        {"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1},
        {"name": "GGG3", "recognition": "GGG", "cut_offset": 0},
        {"name": "CAT3", "recognition": "CAT", "cut_offset": 0},
    ],
    "edits": [
        {"position": 3, "base": "C"},
        {"position": 8, "base": "T"},
    ],
})
check("multi-enzyme report", r["status"] == "solution"
      and r["sequence"] == "GAACTCGGT"
      and r["edit_indices"] == [0, 1] and r["edit_count"] == 2
      and [e["name"] for e in r["enzyme_sites"]] == ["EcoRI", "GGG3", "CAT3"]
      and [e["original_sites"] for e in r["enzyme_sites"]] == [1, 1, 0]
      and [e["remaining_sites"] for e in r["enzyme_sites"]] == [0, 0, 0],
      str(r.get("enzyme_sites", r)))

# 任一目标酶未清零 -> no_solution：只破坏 EcoRI，GGG@6 保留
r = run_erase_sites({
    "sequence": "GAATTCGGG",
    "topology": "linear",
    "enzymes": [
        ECORI,
        {"name": "GGG3", "recognition": "GGG", "cut_offset": 0},
    ],
    "edits": [{"position": 0, "base": "T"}],
})
check("all enzymes must be cleared", r == {"status": "no_solution"}, str(r))

# --- 原本无任何位点：空集合（0 改动）即最优，即使给了候选编辑 ---
r = run_erase_sites({
    "sequence": "ACGTACGT", "topology": "linear", "enzymes": [ECORI],
    "edits": [{"position": 0, "base": "T"}],
})
check("already clean -> zero edits", r["status"] == "solution"
      and r["edit_indices"] == [] and r["edit_count"] == 0
      and r["sequence"] == "ACGTACGT"
      and r["enzyme_sites"][0]["original_sites"] == 0
      and r["enzyme_sites"][0]["remaining_sites"] == 0, str(r))

# --- 无候选编辑且存在位点 -> no_solution ---
r = run_erase_sites({
    "sequence": "AAGAATTCAA", "topology": "linear",
    "enzymes": [ECORI], "edits": [],
})
check("no edits with site -> no_solution", r == {"status": "no_solution"}, str(r))

# --- 环状：编辑后必须重新全环扫描（含跨原点） ---
# 酶识别 GCAC（rc=GTGC）。环 GCGCATGT（len 8）：跨原点位点 GCAC@6
# （@6,7 + @0,1 = G,C,G,C；线性时不存在）。
wrap_enz = {"name": "GC4", "recognition": "GCAC", "cut_offset": 0}
check("circular bait setup one wrap site",
      reference_sites("GCGCATGT", ["GCAC"], "circular")[0] == {6})
# 诱饵编辑 position 1 C->T：跨原点窗口 @6 由 GCAC 翻成 GTGC（反向互补仍命中）。
# 这也验证“同一起点正/反链命中计为一个位点”——翻链后位点依旧存在，并未清除。
bait = run_erase_sites({
    "sequence": "GCGCATGT", "topology": "circular", "enzymes": [wrap_enz],
    "edits": [{"position": 1, "base": "T"}],
})
check("circular strand-flipped wrap site remains -> no_solution",
      bait == {"status": "no_solution"}, str(bait))
# 干净编辑 position 6 G->A：GCGCATAT 环上 GCAC/GTGC 均无
clean = run_erase_sites({
    "sequence": "GCGCATGT", "topology": "circular", "enzymes": [wrap_enz],
    "edits": [
        {"position": 1, "base": "T"},  # 诱饵（下标 0）
        {"position": 6, "base": "A"},  # 干净单编辑（下标 1）
    ],
})
check("circular clean edit chosen over bait",
      clean["status"] == "solution" and clean["edit_indices"] == [1]
      and clean["sequence"] == "GCGCATAT"
      and clean["enzyme_sites"][0]["original_sites"] == 1
      and clean["enzyme_sites"][0]["remaining_sites"] == 0, str(clean))

# 编辑也可能直接补出跨原点位点（原本干净）：
# 环 ACGCATGT（len 8）无 GCAC/GTGC；position 0 A->G 后变 GCGCATGT，
# 跨原点 GCAC@6 出现 —— 空编辑（0 改动）已是最优，该编辑自然不被选中。
r = run_erase_sites({
    "sequence": "ACGCATGT", "topology": "circular", "enzymes": [wrap_enz],
    "edits": [{"position": 0, "base": "G"}],
})
check("circular edit creating wrap site skipped (already clean)",
      r["status"] == "solution" and r["edit_indices"] == []
      and r["sequence"] == "ACGCATGT"
      and r["enzyme_sites"][0]["original_sites"] == 0, str(r))

# --- 长度边界：500 接受，501 拒绝 ---
seq500 = "A" * 494 + "GAATTC"
r = run_erase_sites({
    "sequence": seq500, "topology": "linear", "enzymes": [ECORI],
    "edits": [{"position": 494, "base": "T"}],
})
check("500 accepted and erased", r["status"] == "solution"
      and r["length"] == 500 and r["enzyme_sites"][0]["remaining_sites"] == 0,
      str(r)[:200])
expect_error("501 rejected", {
    "sequence": "A" * 501, "topology": "linear", "enzymes": [ECORI], "edits": [],
}, "上限")

# --- 8 个酶接受，9 个拒绝；10 个编辑接受，11 个拒绝 ---
r8 = run_erase_sites({
    "sequence": "ACGT", "topology": "linear",
    "enzymes": [{"name": f"E{i}", "recognition": "AA", "cut_offset": 0}
                for i in range(8)],
    "edits": [],
})
check("8 enzymes accepted", r8["status"] == "solution"
      and len(r8["enzyme_sites"]) == 8, str(r8)[:200])
expect_error("9 enzymes rejected", {
    "sequence": "ACGT", "topology": "linear",
    "enzymes": [ECORI] * 9, "edits": [],
}, "目标酶数量")
expect_error("0 enzymes rejected", {
    "sequence": "ACGT", "topology": "linear",
    "enzymes": [], "edits": [],
}, "目标酶数量")
# 10 个编辑接受：选一个识别序列从不命中的酶，空编辑即解
r10 = run_erase_sites({
    "sequence": "ACGTACGTAC", "topology": "linear",
    "enzymes": [{"name": "E", "recognition": "TTT", "cut_offset": 0}],
    "edits": [{"position": i,
               "base": "T" if "ACGTACGTAC"[i] != "T" else "A"}
              for i in range(10)],
})
check("10 edits accepted (already clean -> zero chosen)",
      r10["status"] == "solution" and r10["edit_indices"] == []
      and r10["enzyme_sites"][0]["original_sites"] == 0, str(r10)[:200])
expect_error("11 edits rejected", {
    "sequence": "A" * 10, "topology": "linear", "enzymes": [ECORI],
    "edits": [{"position": i, "base": "C"} for i in range(11)],
}, "候选编辑")

# --- 非法输入 ---
expect_error("payload not dict", ["x"], "JSON 对象")
expect_error("missing sequence",
             {"topology": "linear", "enzymes": [ECORI], "edits": []}, "sequence")
expect_error("empty sequence",
             {"sequence": "  ", "topology": "linear",
              "enzymes": [ECORI], "edits": []}, "不能为空")
expect_error("illegal base",
             {"sequence": "ACGXN", "topology": "linear",
              "enzymes": [ECORI], "edits": []}, "非法碱基")
expect_error("bad topology",
             {"sequence": "ACGT", "topology": "circle",
              "enzymes": [ECORI], "edits": []}, "拓扑")
expect_error("enzyme bad offset",
             {"sequence": "ACGT", "topology": "linear",
              "enzymes": [{"name": "E", "recognition": "AA", "cut_offset": 9}],
              "edits": []}, "越界")
expect_error("missing edits",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI]},
             "edits")
expect_error("edits not list",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": {}}, "edits")
expect_error("edit not object",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [5]}, "对象")
expect_error("position missing",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"base": "T"}]}, "position")
expect_error("position not int",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": "1", "base": "T"}]}, "position")
expect_error("position float",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 1.5, "base": "T"}]}, "position")
expect_error("position bool",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": True, "base": "T"}]}, "position")
expect_error("position negative",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": -1, "base": "T"}]}, "越界")
expect_error("position equals length",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 4, "base": "T"}]}, "越界")
expect_error("duplicate position",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 2, "base": "T"},
                        {"position": 2, "base": "G"}]}, "重复")
expect_error("base missing",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 0}]}, "base")
expect_error("bad base letter",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 0, "base": "N"}]}, "base")
expect_error("base not string",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 0, "base": 1}]}, "base")
expect_error("same as original base",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 0, "base": "a"}]}, "必须改变原碱基")
expect_error("same base after strip",
             {"sequence": "ACGT", "topology": "linear", "enzymes": [ECORI],
              "edits": [{"position": 1, "base": " c "}]}, "必须改变原碱基")

# base 大小写与首尾空白归一化
r = run_erase_sites({
    "sequence": "GAATTC", "topology": "linear", "enzymes": [ECORI],
    "edits": [{"position": 0, "base": " t "}],
})
check("base normalized", r["status"] == "solution" and r["sequence"][0] == "T",
      str(r))

# --- 随机不变量：实现结果必须与朴素参考完全一致（含跨原点、正反链、择优规则） ---
random.seed(20260926)
for case_i in range(300):
    n = random.randint(1, 24)
    seq = "".join(random.choice("ACGT") for _ in range(n))
    topology = random.choice(["linear", "circular"])
    enzyme_count = random.randint(1, 3)
    enzymes = []
    recognitions = []
    for e in range(enzyme_count):
        L = random.randint(1, 6)
        rec = "".join(random.choice("ACGT") for _ in range(L))
        recognitions.append(rec)
        enzymes.append({"name": f"E{e}", "recognition": rec, "cut_offset": 0})

    # 随机选若干互不重复位置生成合法候选编辑
    positions = random.sample(range(n), min(n, random.randint(0, 6)))
    edits = [
        {"position": p,
         "base": random.choice([b for b in "ACGT" if b != seq[p]])}
        for p in positions
    ]
    payload = {
        "sequence": seq, "topology": topology,
        "enzymes": enzymes, "edits": edits,
    }
    res = run_erase_sites(payload)
    exp_combo, exp_seq = brute_force_optimal(seq, recognitions, topology, edits)

    if exp_combo is None:
        assert res == {"status": "no_solution"}, (case_i, payload, res, exp_combo)
        continue

    # original 计数与参考一致（可行解响应中才有 enzyme_sites）
    ref_original = reference_sites(seq, recognitions, topology)
    assert [e["original_sites"] for e in res["enzyme_sites"]] \
        == [len(s) for s in ref_original], (case_i, payload, res)
    assert res["status"] == "solution", (case_i, payload, res)
    assert res["edit_indices"] == exp_combo, (case_i, payload, res, exp_combo)
    assert res["sequence"] == exp_seq, (case_i, payload, res, exp_seq)
    assert len(res["sequence"]) == n
    ref_remaining = reference_sites(exp_seq, recognitions, topology)
    assert [e["remaining_sites"] for e in res["enzyme_sites"]] \
        == [len(s) for s in ref_remaining]
    assert all(len(s) == 0 for s in ref_remaining)
    # 编辑只发生在所选位置，其余碱基不变
    changed = {edits[i]["position"] for i in exp_combo}
    for p in range(n):
        if p not in changed:
            assert res["sequence"][p] == seq[p]
print("  PASS randomized equivalence with brute force (300 cases)")

print(f"\n{PASS + 1} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
