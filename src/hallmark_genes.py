"""
Hallmark gene sets for the 9 hallmarks of aging (López-Otín et al., 2013).
Curated from GenAge, MSigDB aging pathways, and primary literature.
"""

HALLMARK_GENES = {
    "Genomic Instability": [
        "TP53", "ATM", "ATR", "BRCA1", "BRCA2", "RAD51", "XRCC1", "XRCC5",
        "XRCC6", "PARP1", "MLH1", "MSH2", "MSH6", "ERCC1", "ERCC2", "XPC",
        "OGG1", "APEX1", "LIG3", "POLB", "BLM", "WRN", "RECQL4", "LMNA",
        "BUB1B", "CHEK1", "CHEK2", "MDM2", "CDKN1A", "H2AFX", "NBN",
        "MRE11", "RAD50", "RPA1", "PCNA", "RFC1"
    ],
    "Telomere Attrition": [
        "TERT", "TERC", "DKC1", "TINF2", "ACD", "POT1", "TERF1", "TERF2",
        "TERF2IP", "RTEL1", "CTC1", "STN1", "TEN1", "WRAP53", "NHP2",
        "NOP10", "GAR1", "NAF1", "PARN", "OBFC1"
    ],
    "Epigenetic Alterations": [
        "SIRT1", "SIRT3", "SIRT6", "SIRT7", "DNMT1", "DNMT3A", "DNMT3B",
        "TET1", "TET2", "TET3", "HDAC1", "HDAC2", "HDAC3", "KAT2A",
        "KAT2B", "EP300", "CREBBP", "EZH2", "SUV39H1", "KDM1A", "KDM6A",
        "KDM6B", "CBX5", "BMI1", "RING1", "SETDB1", "KMT2A", "KMT2D",
        "DOT1L", "PRMT1"
    ],
    "Loss of Proteostasis": [
        "HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1", "HSPB1",
        "HSF1", "HSF2", "BAG3", "STUB1", "UBB", "UBC", "USP14", "PSMD11",
        "PSMA7", "PSMB5", "SQSTM1", "MAP1LC3B", "BECN1", "ATG5", "ATG7",
        "ATG12", "LAMP2", "CTSD", "CTSB", "TFEB", "ULK1", "VCP", "UBQLN1",
        "UBQLN2"
    ],
    "Deregulated Nutrient Sensing": [
        "IGF1", "IGF1R", "INSR", "IRS1", "IRS2", "AKT1", "AKT2", "MTOR",
        "RPTOR", "RICTOR", "RPS6KB1", "EIF4EBP1", "FOXO1", "FOXO3",
        "FOXO4", "PRKAA1", "PRKAA2", "STK11", "TSC1", "TSC2", "RHEB",
        "GHR", "PTEN", "PIK3CA", "PIK3CB", "PPARGC1A", "PPARGC1B",
        "ADIPOQ", "LEP", "SLC2A4"
    ],
    "Mitochondrial Dysfunction": [
        "PPARGC1A", "PPARGC1B", "TFAM", "NRF1", "POLG", "POLG2", "TWNK",
        "MT-ND1", "MT-ND4", "MT-ND5", "MT-CO1", "MT-CO2", "MT-CO3",
        "MT-ATP6", "MT-ATP8", "MT-CYB", "SOD2", "SOD1", "CAT", "GPX1",
        "GPX4", "PRDX3", "PARK2", "PINK1", "MFN1", "MFN2", "DNM1L",
        "FIS1", "OPA1", "SIRT3", "UCP2", "CYCS"
    ],
    "Cellular Senescence": [
        "CDKN2A", "CDKN2B", "CDKN1A", "TP53", "RB1", "MDM2", "CCND1",
        "CDK4", "CDK6", "E2F1", "SERPINE1", "IL6", "IL8", "CXCL1",
        "CXCL2", "CCL2", "MMP1", "MMP3", "MMP9", "IGFBP3", "IGFBP5",
        "IGFBP7", "LMNB1", "GLB1", "HMGA1", "HMGA2", "H2AFX", "TGFB1",
        "BMP2", "GDF15"
    ],
    "Stem Cell Exhaustion": [
        "NANOG", "POU5F1", "SOX2", "KLF4", "MYC", "BMI1", "EZH2",
        "TERT", "CDKN2A", "CDKN1A", "TP53", "NOTCH1", "NOTCH2", "WNT3A",
        "CTNNB1", "LGR5", "PROM1", "THY1", "ALDH1A1", "FGF2", "FGFR1",
        "HGF", "KIT", "SHH", "BMP4", "LIF", "STAT3", "JAK2", "CDC42",
        "FOXO3"
    ],
    "Altered Intercellular Communication": [
        "NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA", "IKBKB", "NLRP3",
        "PYCARD", "CASP1", "IL1B", "IL6", "TNF", "IL18", "IL10", "TGFB1",
        "IFNG", "CCL2", "CCL5", "CXCL12", "CRP", "GNRH1", "GNRHR",
        "TLR2", "TLR4", "MYD88", "TICAM1", "HMGB1", "AGER", "SERPINE1",
        "VEGFA"
    ],
}

# Hierarchical classification
HALLMARK_CATEGORIES = {
    "Primary": [
        "Genomic Instability",
        "Telomere Attrition",
        "Epigenetic Alterations",
        "Loss of Proteostasis",
    ],
    "Antagonistic": [
        "Deregulated Nutrient Sensing",
        "Mitochondrial Dysfunction",
        "Cellular Senescence",
    ],
    "Integrative": [
        "Stem Cell Exhaustion",
        "Altered Intercellular Communication",
    ],
}

HALLMARK_SHORT = {
    "Genomic Instability": "GI",
    "Telomere Attrition": "TA",
    "Epigenetic Alterations": "EA",
    "Loss of Proteostasis": "LP",
    "Deregulated Nutrient Sensing": "DNS",
    "Mitochondrial Dysfunction": "MD",
    "Cellular Senescence": "CS",
    "Stem Cell Exhaustion": "SCE",
    "Altered Intercellular Communication": "AIC",
}


def get_all_genes():
    """Return a deduplicated sorted list of all hallmark genes."""
    genes = set()
    for gene_list in HALLMARK_GENES.values():
        genes.update(gene_list)
    return sorted(genes)


def get_gene_hallmark_map():
    """Return a dict mapping each gene to its hallmark(s)."""
    mapping = {}
    for hallmark, genes in HALLMARK_GENES.items():
        for gene in genes:
            mapping.setdefault(gene, []).append(hallmark)
    return mapping


def get_shared_genes():
    """Return genes that appear in multiple hallmarks (cross-talk nodes)."""
    gene_map = get_gene_hallmark_map()
    return {g: h for g, h in gene_map.items() if len(h) > 1}
