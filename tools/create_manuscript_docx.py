from __future__ import annotations

import zipfile
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "TA4_work_hardening_manuscript_draft.docx"


CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.settings+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"""

RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"""

WORD_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>
"""

APP = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Microsoft Word</Application>
  <DocSecurity>0</DocSecurity>
  <ScaleCrop>false</ScaleCrop>
  <Company></Company>
  <LinksUpToDate>false</LinksUpToDate>
  <SharedDoc>false</SharedDoc>
  <HyperlinksChanged>false</HyperlinksChanged>
  <AppVersion>16.0000</AppVersion>
</Properties>
"""

CORE = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>TA4 Work Hardening Manuscript Draft</dc:title>
  <dc:creator>Codex</dc:creator>
  <cp:lastModifiedBy>Codex</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">2026-05-21T00:00:00Z</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">2026-05-21T00:00:00Z</dcterms:modified>
</cp:coreProperties>
"""

SETTINGS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:settings xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:zoom w:percent="100"/>
  <w:defaultTabStop w:val="420"/>
  <w:characterSpacingControl w:val="doNotCompress"/>
  <w:compat/>
</w:settings>
"""

STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:docDefaults>
    <w:rPrDefault>
      <w:rPr>
        <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="宋体"/>
        <w:sz w:val="24"/>
        <w:szCs w:val="24"/>
        <w:lang w:val="en-US" w:eastAsia="zh-CN"/>
      </w:rPr>
    </w:rPrDefault>
    <w:pPrDefault>
      <w:pPr>
        <w:spacing w:after="120" w:line="360" w:lineRule="auto"/>
      </w:pPr>
    </w:pPrDefault>
  </w:docDefaults>
  <w:style w:type="paragraph" w:default="1" w:styleId="Normal">
    <w:name w:val="Normal"/>
    <w:qFormat/>
    <w:pPr>
      <w:spacing w:after="120" w:line="360" w:lineRule="auto"/>
      <w:jc w:val="both"/>
    </w:pPr>
    <w:rPr>
      <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="宋体"/>
      <w:sz w:val="24"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Title">
    <w:name w:val="Title"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr>
      <w:spacing w:before="240" w:after="240" w:line="360" w:lineRule="auto"/>
      <w:jc w:val="center"/>
    </w:pPr>
    <w:rPr>
      <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="黑体"/>
      <w:b/>
      <w:sz w:val="32"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Subtitle">
    <w:name w:val="Subtitle"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr>
      <w:spacing w:after="240" w:line="320" w:lineRule="auto"/>
      <w:jc w:val="center"/>
    </w:pPr>
    <w:rPr>
      <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="宋体"/>
      <w:i/>
      <w:sz w:val="22"/>
      <w:color w:val="555555"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading1">
    <w:name w:val="heading 1"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr>
      <w:keepNext/>
      <w:spacing w:before="360" w:after="160" w:line="360" w:lineRule="auto"/>
      <w:outlineLvl w:val="0"/>
    </w:pPr>
    <w:rPr>
      <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="黑体"/>
      <w:b/>
      <w:sz w:val="28"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Heading2">
    <w:name w:val="heading 2"/>
    <w:basedOn w:val="Normal"/>
    <w:next w:val="Normal"/>
    <w:qFormat/>
    <w:pPr>
      <w:keepNext/>
      <w:spacing w:before="240" w:after="120" w:line="320" w:lineRule="auto"/>
      <w:outlineLvl w:val="1"/>
    </w:pPr>
    <w:rPr>
      <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="黑体"/>
      <w:b/>
      <w:sz w:val="24"/>
    </w:rPr>
  </w:style>
  <w:style w:type="paragraph" w:styleId="Placeholder">
    <w:name w:val="Placeholder"/>
    <w:basedOn w:val="Normal"/>
    <w:pPr>
      <w:spacing w:after="160" w:line="360" w:lineRule="auto"/>
    </w:pPr>
    <w:rPr>
      <w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="宋体"/>
      <w:color w:val="666666"/>
      <w:sz w:val="23"/>
    </w:rPr>
  </w:style>
</w:styles>
"""


def paragraph(text: str = "", style: str = "Normal", page_break_before: bool = False) -> str:
    ppr = f"<w:pStyle w:val=\"{style}\"/>"
    if page_break_before:
        ppr += "<w:pageBreakBefore/>"
    escaped = escape(text)
    return (
        "<w:p>"
        f"<w:pPr>{ppr}</w:pPr>"
        f"<w:r><w:t xml:space=\"preserve\">{escaped}</w:t></w:r>"
        "</w:p>"
    )


def blank() -> str:
    return "<w:p><w:pPr><w:pStyle w:val=\"Normal\"/></w:pPr></w:p>"


def make_document() -> str:
    paragraphs: list[str] = []
    paragraphs.append(paragraph("旋锻冷变形对 TA4 商业纯钛加工硬化行为与组织机制的影响", "Title"))
    paragraphs.append(paragraph("参考《7075 高强铝合金构件冷成形强化机制研究》的期刊论文结构", "Subtitle"))
    paragraphs.append(paragraph("作者：__________    单位：__________    日期：__________", "Placeholder"))
    paragraphs.append(blank())
    paragraphs.append(paragraph("摘要", "Heading1"))
    paragraphs.append(paragraph("目的  【待写】针对 TA4 商业纯钛冷变形强化过程中强度提升、塑性变化及加工硬化机制尚需进一步阐明的问题，说明本文研究目的。", "Placeholder"))
    paragraphs.append(paragraph("方法  【待写】说明采用不同旋锻变形量样品，结合拉伸试验、真实应力-应变分析、加工硬化率计算、金相组织观察和 EBSD 表征开展研究。", "Placeholder"))
    paragraphs.append(paragraph("结果  【待写】概述屈服强度、抗拉强度、延伸率和加工硬化行为随冷变形量的变化，并说明主要金相/EBSD 观察结果。", "Placeholder"))
    paragraphs.append(paragraph("结论  【待写】归纳 TA4 商业纯钛冷变形强化的主要机制，注意将位错积累、晶界/亚晶界作用、取向或织构演化等结论限定在数据支持范围内。", "Placeholder"))
    paragraphs.append(paragraph("关键词：TA4；商业纯钛；旋锻；冷变形；加工硬化；EBSD", "Placeholder"))
    paragraphs.append(blank())
    paragraphs.append(paragraph("Work Hardening Behavior and Microstructural Mechanism of Rotary-swaged TA4 Commercially Pure Titanium", "Heading1"))
    paragraphs.append(paragraph("ABSTRACT: 【To be written after the Chinese abstract is finalized. Keep the same purpose-method-results-conclusion logic as the Chinese abstract.】", "Placeholder"))
    paragraphs.append(paragraph("KEY WORDS: TA4; commercially pure titanium; rotary swaging; cold deformation; work hardening; EBSD", "Placeholder"))

    sections = [
        ("引言", [
            ("", "【待写】参考目标文章的写法，本部分不单独分成多个小节。建议按“材料应用背景和性能需求 -> 冷变形/旋锻强化研究现状 -> 现有研究不足 -> 本文研究对象、方法和目的”的顺序展开。"),
            ("", "【待写】可写明：本文以 TA4 商业纯钛为研究对象，分析不同旋锻冷变形量下拉伸性能、加工硬化行为和显微组织演化之间的关系，以阐明冷变形诱导强度提升及塑性变化的组织机制。"),
        ]),
        ("1 试验", [
            ("", "选取 12 根初始直径为 7 mm 的 TA4 商业纯钛棒材为研究对象。其中 2 根棒材保留原始状态，记为 TA4-M；其余棒材经室温旋锻冷变形加工至不同直径，记为 TA4-Y-d，其中 d 表示旋锻后的名义直径，单位为 mm。依据实验室拉伸数据中的直径分组，旋锻后样品直径分别为 6.5、6.02、5.6、5.25 和 5.0 mm，各直径组包含 2 根平行样品。因此，旋锻态样品分别记为 TA4-Y-6.5、TA4-Y-6.02、TA4-Y-5.6、TA4-Y-5.25 和 TA4-Y-5.0。"),
            ("", "冷变形量按棒材截面积缩减率计算，如式（1）所示："),
            ("", "ε = (A0 - A1) / A0 × 100% = (D0^2 - D1^2) / D0^2 × 100%    （1）"),
            ("", "式中，A0 和 A1 分别为旋锻前、后的横截面积，D0 和 D1 分别为旋锻前、后的棒材直径。以 D0 = 7 mm 作为参考直径计算，各样品组的冷变形量分别为 0、13.78%、26.04%、36.00%、43.75% 和 48.98%。"),
            ("", "对不同状态的棒材进行室温单轴拉伸试验，获得工程应力-应变曲线，并统计 0.2% 偏移屈服强度、抗拉强度、断后伸长率和断面收缩率。每个直径组设置 2 个平行试样，结果以平均值表示；对于弹性段不清晰、0.2% 偏移屈服点无法可靠确定或曲线异常的结果，在后续数据分析中保留可靠性标记，并结合原始曲线进行判别。"),
            ("", "为分析旋锻冷变形对 TA4 商业纯钛组织演化的影响，分别选取 TA4-M 及不同 TA4-Y-d 样品进行金相观察和 EBSD 表征。金相样品经切取、镶嵌、磨抛和腐蚀后，用于观察不同变形量下晶粒形貌、变形流线和组织均匀性。EBSD 表征用于获得不同状态样品的取向分布、晶界特征、局部取向差及织构演化信息，并与拉伸性能和加工硬化行为相结合，分析旋锻冷变形诱导的强化机制。"),
        ]),
        ("2 结果与分析", [
            ("2.1 不同冷变形量下 TA4 的力学性能", "【待写】对应参考文章的“板料力学性能”小节。展示不同直径/冷变形量样品的屈服强度、抗拉强度和延伸率。可先给出应力-应变曲线，再给出性能统计图或表。"),
            ("", "【待写】对于偏离整体变化趋势的结果，使用“非单调变化”“局部差异”“离散性增加”等学术表述，并结合原始曲线可靠性、试样状态和组织证据分析。"),
            ("2.2 真实应力-应变曲线与加工硬化行为", "【待写】对应参考文章的“成形构件力学性能”小节，但根据本课题改为加工硬化行为。展示真实应力-真实应变曲线、加工硬化率 theta = d sigma / d epsilon，以及必要时的归一化加工硬化率或 Hollomon 指数。"),
            ("", "【待写】重点分析冷变形量提高后材料初始强度、均匀变形能力和后续加工硬化空间之间的关系。"),
            ("2.3 TA4 商业纯钛冷变形强化机制", "【待写】对应参考文章的“冷成形强化机制”小节，是全文机制讨论的核心。结合金相和 EBSD 结果，讨论位错积累、晶界/亚晶界阻碍、取向梯度或织构演化对强度提升和塑性变化的影响。"),
            ("", "【待写】若 EBSD 支持，可讨论 LAGB/HAGB、KAM、取向差分布和织构变化。若无直接证据，不应定量声称 GND 密度或孪生机制。"),
        ]),
        ("3 结论", [
            ("", "1）【待写】概括旋锻冷变形量对 TA4 商业纯钛屈服强度、抗拉强度和延伸率的影响。"),
            ("", "2）【待写】概括不同冷变形量下真实应力-应变曲线和加工硬化率的主要变化特征。"),
            ("", "3）【待写】概括金相和 EBSD 所支持的冷变形强化机制，注意将解释限定在位错积累、晶界/亚晶界、取向/织构演化等已有证据范围内。"),
            ("", "4）【待写】如存在偏离整体趋势的结果，简要说明可能原因及仍需进一步验证的组织指标。"),
        ]),
        ("参考文献：", [
            ("", "【待补】按目标期刊格式整理。正文中每一处机制解释应对应实验数据或文献依据。"),
        ]),
    ]

    for title, children in sections:
        paragraphs.append(paragraph(title, "Heading1"))
        for subtitle, text in children:
            if subtitle:
                paragraphs.append(paragraph(subtitle, "Heading2"))
            paragraphs.append(paragraph(text, "Placeholder"))
            paragraphs.append(blank())

    sect_pr = (
        "<w:sectPr>"
        "<w:pgSz w:w=\"12240\" w:h=\"15840\"/>"
        "<w:pgMar w:top=\"1440\" w:right=\"1440\" w:bottom=\"1440\" w:left=\"1440\" w:header=\"720\" w:footer=\"720\" w:gutter=\"0\"/>"
        "<w:cols w:space=\"720\"/>"
        "<w:docGrid w:linePitch=\"360\"/>"
        "</w:sectPr>"
    )
    body = "".join(paragraphs) + sect_pr
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
        'xmlns:xml="http://www.w3.org/XML/1998/namespace">'
        f"<w:body>{body}</w:body>"
        "</w:document>"
    )


def main() -> None:
    with zipfile.ZipFile(OUT, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", CONTENT_TYPES)
        zf.writestr("_rels/.rels", RELS)
        zf.writestr("docProps/app.xml", APP)
        zf.writestr("docProps/core.xml", CORE)
        zf.writestr("word/_rels/document.xml.rels", WORD_RELS)
        zf.writestr("word/document.xml", make_document())
        zf.writestr("word/styles.xml", STYLES)
        zf.writestr("word/settings.xml", SETTINGS)
    print(OUT)


if __name__ == "__main__":
    main()
