import pandas as pd
import numpy as np
from cyvcf2 import VCF
from collections import defaultdict

# Annotating P, P/LP and VUS Clinvar variants with known rsIDs in GI
clinvar_gi_id_ref_alt = pd.read_csv("rsid_known/id_ref_alt.txt", sep="\t", names=["SNP", "Ref", "Alt"])
clinvar_gi_id_ref_alt["SNP"] = clinvar_gi_id_ref_alt["SNP"].str.split(";")
clinvar_gi_id_ref_alt = clinvar_gi_id_ref_alt.explode("SNP").reset_index(drop=True)

clinvar_gi_freq = pd.read_csv("rsid_known/clinvar_p_lp_vus_gi_using_rsids.frq", sep=r"\s+", usecols=["SNP", "A1", "A2", "MAF"])
clinvar_gi_freq["SNP"] = clinvar_gi_freq["SNP"].str.split(";")
clinvar_gi_freq = clinvar_gi_freq.explode("SNP").reset_index(drop=True)

clinvar_gi_allele_count = pd.read_csv("rsid_known/allele_count.frq.counts", sep=r"\s+", usecols=["SNP", "C1", "C2"])
clinvar_gi_allele_count["SNP"] = clinvar_gi_allele_count["SNP"].str.split(";")
clinvar_gi_allele_count = clinvar_gi_allele_count.explode("SNP").reset_index(drop=True)

clinvar_gi_genotype_count = pd.read_csv("rsid_known/genotype_count.frqx", sep="\t", usecols=["SNP", "C(HOM A1)", "C(HET)", "C(HOM A2)"])
clinvar_gi_genotype_count["SNP"] = clinvar_gi_genotype_count["SNP"].str.split(";")
clinvar_gi_genotype_count = clinvar_gi_genotype_count.explode("SNP").reset_index(drop=True)

clinvar_p_lp_vus_with_rsid = pd.read_csv("rsid_known/clinvar_with_rsid.txt", sep="\t")

acmg_genes = pd.read_excel("ACMG_gene_list.xlsx", header=2)
acmg_genes_filtered = acmg_genes[acmg_genes["Gene MIM"].notna()]
genes = acmg_genes_filtered["Gene"].drop_duplicates()
acmg_gene_list = genes.to_list()

clinvar_ref_alt_freq = clinvar_gi_id_ref_alt.merge(clinvar_gi_freq, on="SNP").merge(clinvar_gi_allele_count, on="SNP").merge(clinvar_gi_genotype_count, on="SNP")
clinvar_with_rsid_merge = pd.merge(clinvar_ref_alt_freq, clinvar_p_lp_vus_with_rsid, left_on="SNP", right_on="rsRS# (dbSNP)", how="inner")

final_clinvar_rsid = clinvar_with_rsid_merge.rename(columns={"SNP":"variant_ID", "A1":"GI_Minor_Allele", "A2":"GI_Major_Allele", "Ref":"GI_Ref", "Alt":"GI_Alt", 
                                                             "MAF":"GI_MAF", "C1":"GI_Minor_Allele_Count", "C2":"GI_Major_Allele_Count", "C(HOM A1)":"GI_HOM_Minor_Allele_Count",
                                                             "C(HET)":"GI_HET", "C(HOM A2)": "GI_HOM_Major_Allele_Count", "ReferenceAlleleVCF":"clinvar_ref",
                                                             "AlternateAlleleVCF":"clinvar_alt"})

final_clinvar_rsid["ACMG_Gene"] = final_clinvar_rsid["GeneSymbol"].isin(acmg_gene_list).map({True:"Yes", False:"No"})
def rev_comp(allele):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp.get(base, base) for base in allele)

final_clinvar_rsid["GI_Clinvar_Match"] = np.where(
    (((final_clinvar_rsid["GI_Ref"] == final_clinvar_rsid["clinvar_ref"]) & 
    (final_clinvar_rsid["GI_Alt"] == final_clinvar_rsid["clinvar_alt"])) |
    ((final_clinvar_rsid["GI_Ref"].apply(rev_comp) == final_clinvar_rsid["clinvar_ref"]) &
     (final_clinvar_rsid["GI_Alt"].apply(rev_comp) == final_clinvar_rsid["clinvar_alt"]))),
     "Yes",
     "No"
)
final_clinvar_rsid = final_clinvar_rsid.fillna("-")
final_clinvar_rsid["Chromosome"] = pd.to_numeric(final_clinvar_rsid["Chromosome"], errors="coerce")

final_clinvar_rsid = final_clinvar_rsid.drop_duplicates()
final_clinvar_rsid = final_clinvar_rsid[["variant_ID", "Chromosome", "Start", "Stop", "Name", "Type", "GeneSymbol", "ACMG_Gene", 
                                         "clinvar_ref", "clinvar_alt", "GI_Ref", "GI_Alt", "GI_Major_Allele", "GI_Minor_Allele", 
                                         "GI_Clinvar_Match", "Origin", "ClinicalSignificance", "PhenotypeList", "GI_MAF", "GI_Major_Allele_Count",
                                         "GI_Minor_Allele_Count", "GI_HOM_Major_Allele_Count", "GI_HET", "GI_HOM_Minor_Allele_Count"]]

final_clinvar_rsid.to_excel("rsid_known/clinvar_gi_with_rsid.xlsx", index=False)

# Annotating P, P/LP and VUS Clinvar variants with unknown rsIDs in GI
clinvar_p_lp_vus_no_rsid = pd.read_csv("rsid_unknown_clinvar/clinvar_no_rsid.txt", sep="\t", names=['Type', 'Name', 'GeneSymbol', 
                                                                                                    'ClinicalSignificance', 
                                                                                                    'rsRS# (dbSNP)', 'PhenotypeList', 'Origin', 
                                                                                                    'Assembly', 'Chromosome', 
                                                                                                    'Start', 'Stop', 'ReferenceAlleleVCF', 'AlternateAlleleVCF'])
snvs_only = clinvar_p_lp_vus_no_rsid[clinvar_p_lp_vus_no_rsid["Type"] == "single nucleotide variant"]
not_snvs = clinvar_p_lp_vus_no_rsid[
    (clinvar_p_lp_vus_no_rsid["Type"] != "single nucleotide variant") &
    (clinvar_p_lp_vus_no_rsid["ReferenceAlleleVCF"] != "na") &
    (clinvar_p_lp_vus_no_rsid["AlternateAlleleVCF"] != "na")
]

no_rsid_ref_alt_gi = pd.read_csv("rsid_unknown_clinvar/snp_ref_alt.txt", sep="\t", names=['Chromosome', 'Position', 'variant_ID', 
                                                                                          'GI_Ref', 'GI_Alt', 'GI_MAF'])
no_rsid_ref_alt_gi["variant_ID"] = no_rsid_ref_alt_gi["variant_ID"].str.split(";")
no_rsid_ref_alt_gi = no_rsid_ref_alt_gi.explode("variant_ID").reset_index(drop=True)
no_rsid_ref_alt_gi['Chromosome'] = no_rsid_ref_alt_gi['Chromosome'].str.replace('chr', '', regex=False)
no_rsid_ref_alt_gi['Chromosome'] = no_rsid_ref_alt_gi['Chromosome'].astype(str)

indel_ref_alt_gi = pd.read_csv("rsid_unknown_clinvar/indel_ref_alt.txt", sep="\t", names=['Chromosome', 'Position', 'variant_ID', 
                                                                                          'GI_Ref', 'GI_Alt', 'GI_MAF'])
indel_ref_alt_gi["variant_ID"] = indel_ref_alt_gi["variant_ID"].str.split(";")
indel_ref_alt_gi = indel_ref_alt_gi.explode("variant_ID").reset_index(drop=True)
indel_ref_alt_gi['Chromosome'] = indel_ref_alt_gi['Chromosome'].str.replace('chr', '', regex=False)
indel_ref_alt_gi['Chromosome'] = indel_ref_alt_gi['Chromosome'].astype(str)

no_rsid_bim = pd.read_csv("rsid_unknown_clinvar/snp_no_id.bim", sep="\t", names=['Chromosome', 'ID_bim', 'cM', 'Position', 
                                                                                 'GI_Minor_Allele', 'GI_Major_Allele'])
no_rsid_bim["ID_bim"] = no_rsid_bim["ID_bim"].str.split(";")
no_rsid_bim = no_rsid_bim.explode("ID_bim").reset_index(drop=True)
no_rsid_bim = no_rsid_bim[['GI_Minor_Allele', 'GI_Major_Allele']]

indel_bim = pd.read_csv("rsid_unknown_clinvar/indel_no_id.bim", sep="\t", names=['Chromosome', 'ID_bim', 'cM', 'Position', 
                                                                                 'GI_Minor_Allele', 'GI_Major_Allele'])
indel_bim["ID_bim"] = indel_bim["ID_bim"].str.split(";")
indel_bim = indel_bim.explode("ID_bim").reset_index(drop=True)
indel_bim = indel_bim[['GI_Minor_Allele', 'GI_Major_Allele']]

snv_allele_count = pd.read_csv("rsid_unknown_clinvar/allele_count.frq.counts", sep=r"\s+", usecols=["SNP","A1", "A2", "C1", "C2"])
snv_allele_count["SNP"] = snv_allele_count["SNP"].str.split(";")
snv_allele_count = snv_allele_count.explode("SNP").reset_index(drop=True)
snv_allele_count = snv_allele_count[["C1", "C2"]]

indel_allele_count = pd.read_csv("rsid_unknown_clinvar/indel_allele_count.frq.counts", sep=r"\s+", usecols=["SNP","A1", "A2", "C1", "C2"])
indel_allele_count["SNP"] = indel_allele_count["SNP"].str.split(";")
indel_allele_count = indel_allele_count.explode("SNP").reset_index(drop=True)
indel_allele_count = indel_allele_count[["C1", "C2"]]

snv_genotype_count = pd.read_csv("rsid_unknown_clinvar/genotype_count.frqx", sep="\t", usecols=["SNP", "C(HOM A1)", "C(HET)", "C(HOM A2)"])
snv_genotype_count["SNP"] = snv_genotype_count["SNP"].str.split(";")
snv_genotype_count = snv_genotype_count.explode("SNP").reset_index(drop=True)
snv_genotype_count = snv_genotype_count[["C(HOM A1)", "C(HET)", "C(HOM A2)"]]

indel_genotype_count = pd.read_csv("rsid_unknown_clinvar/indel_genotype_count.frqx", sep="\t", usecols=["SNP", "C(HOM A1)", "C(HET)", "C(HOM A2)"])
indel_genotype_count["SNP"] = indel_genotype_count["SNP"].str.split(";")
indel_genotype_count = indel_genotype_count.explode("SNP").reset_index(drop=True)
indel_genotype_count = indel_genotype_count[["C(HOM A1)", "C(HET)", "C(HOM A2)"]]

merge_bim_ref_alt_snv = pd.concat([no_rsid_ref_alt_gi, no_rsid_bim, snv_allele_count, snv_genotype_count], axis=1)
merge_bim_ref_alt_indel = pd.concat([indel_ref_alt_gi, indel_bim, indel_allele_count, indel_genotype_count], axis=1)

final_clinvar_no_rsid_snv = pd.merge(merge_bim_ref_alt_snv, snvs_only,
                                     left_on=['Chromosome', 'Position', 'Position'],
                                     right_on=['Chromosome', 'Start', 'Stop'],
                                     how='inner')

final_clinvar_no_rsid_indel = pd.merge(merge_bim_ref_alt_indel, not_snvs,
                                       left_on=['Chromosome', 'Position'],
                                       right_on=['Chromosome', 'Start'],
                                       how='inner')

final_clinvar_no_rsid_snv["ACMG_Gene"] = final_clinvar_no_rsid_snv["GeneSymbol"].isin(acmg_gene_list).map({True:"Yes", False:"No"})
final_clinvar_no_rsid_indel["ACMG_Gene"] = final_clinvar_no_rsid_indel["GeneSymbol"].isin(acmg_gene_list).map({True:"Yes", False:"No"})

final_clinvar_no_rsid_snv = final_clinvar_no_rsid_snv.rename(columns={"C1":"GI_Minor_Allele_Count", "C2":"GI_Major_Allele_Count", "C(HOM A1)":"GI_HOM_Minor_Allele_Count",
                                                              "C(HET)":"GI_HET", "C(HOM A2)":"GI_HOM_Major_Allele_Count", "ReferenceAlleleVCF":"clinvar_ref",
                                                              "AlternateAlleleVCF":"clinvar_alt"})

final_clinvar_no_rsid_indel = final_clinvar_no_rsid_indel.rename(columns={"C1":"GI_Minor_Allele_Count", "C2":"GI_Major_Allele_Count", "C(HOM A1)":"GI_HOM_Minor_Allele_Count",
                                                              "C(HET)":"GI_HET", "C(HOM A2)":"GI_HOM_Major_Allele_Count", "ReferenceAlleleVCF":"clinvar_ref",
                                                              "AlternateAlleleVCF":"clinvar_alt"})

final_clinvar_no_rsid_snv["GI_Clinvar_Match"] = np.where(
    (((final_clinvar_no_rsid_snv["GI_Ref"] == final_clinvar_no_rsid_snv["clinvar_ref"]) & 
    (final_clinvar_no_rsid_snv["GI_Alt"] == final_clinvar_no_rsid_snv["clinvar_alt"])) |
    ((final_clinvar_no_rsid_snv["GI_Ref"].apply(rev_comp) == final_clinvar_no_rsid_snv["clinvar_ref"]) &
    (final_clinvar_no_rsid_snv["GI_Alt"].apply(rev_comp) == final_clinvar_no_rsid_snv["clinvar_alt"]))),
    "Yes",
    "No"
)

final_clinvar_no_rsid_indel["GI_Clinvar_Match"] = np.where(
    (((final_clinvar_no_rsid_indel["GI_Ref"] == final_clinvar_no_rsid_indel["clinvar_ref"]) & 
    (final_clinvar_no_rsid_indel["GI_Alt"] == final_clinvar_no_rsid_indel["clinvar_alt"])) |
    ((final_clinvar_no_rsid_indel["GI_Ref"].apply(rev_comp) == final_clinvar_no_rsid_indel["clinvar_ref"]) &
    (final_clinvar_no_rsid_indel["GI_Alt"].apply(rev_comp) == final_clinvar_no_rsid_indel["clinvar_alt"]))),
    "Yes",
    "No"
)

final_clinvar_no_rsid_snv = final_clinvar_no_rsid_snv.fillna("-")
final_clinvar_no_rsid_snv["Chromosome"] = pd.to_numeric(final_clinvar_no_rsid_snv["Chromosome"], errors="coerce")

final_clinvar_no_rsid_indel = final_clinvar_no_rsid_indel.fillna("-")
final_clinvar_no_rsid_indel["Chromosome"] = pd.to_numeric(final_clinvar_no_rsid_indel["Chromosome"], errors="coerce")

final_clinvar_no_rsid_snv = final_clinvar_no_rsid_snv.drop_duplicates()
final_clinvar_no_rsid_snv = final_clinvar_no_rsid_snv[["variant_ID", "Chromosome", "Start", "Stop", "Name", "Type", "GeneSymbol", 
                                               "ACMG_Gene", "clinvar_ref", "clinvar_alt", "GI_Ref", "GI_Alt", "GI_Major_Allele", 
                                               "GI_Minor_Allele","GI_Clinvar_Match", "Origin", "ClinicalSignificance", "PhenotypeList", "GI_MAF",
                                               "GI_Major_Allele_Count", "GI_Minor_Allele_Count", "GI_HOM_Major_Allele_Count", "GI_HET", "GI_HOM_Minor_Allele_Count"]]

final_clinvar_no_rsid_indel = final_clinvar_no_rsid_indel.drop_duplicates()
final_clinvar_no_rsid_indel = final_clinvar_no_rsid_indel[["variant_ID", "Chromosome", "Start", "Stop", "Name", "Type", "GeneSymbol", 
                                               "ACMG_Gene", "clinvar_ref", "clinvar_alt", "GI_Ref", "GI_Alt", "GI_Major_Allele", 
                                               "GI_Minor_Allele","GI_Clinvar_Match", "Origin", "ClinicalSignificance", "PhenotypeList", "GI_MAF",
                                               "GI_Major_Allele_Count", "GI_Minor_Allele_Count", "GI_HOM_Major_Allele_Count", "GI_HET", "GI_HOM_Minor_Allele_Count"]]

final_clinvar_no_rsid_snv.to_excel("rsid_unknown_clinvar/clinvar_gi_with_no_rsid_snv.xlsx", index=False)
final_clinvar_no_rsid_indel.to_excel("rsid_unknown_clinvar/clinvar_gi_with_no_rsid_indel.xlsx", index=False)

final_clinvar_gi = pd.concat([final_clinvar_rsid, final_clinvar_no_rsid_snv, final_clinvar_no_rsid_indel], axis=0, ignore_index=True)
final_clinvar_gi = final_clinvar_gi.drop_duplicates()
final_clinvar_gi.to_excel("final_clinvar_gi.xlsx", index=False)

# Per-sample genotype class of each Clinvar P, P/LP and VUS variants in GI
vcf = VCF("var_pos_gi_unrequired_ids_removed_unrequired_pos_removed.vcf.gz")
samples = vcf.samples

rows = []

for var in vcf:
    for sample, gt in zip(samples, var.genotypes):
        # gt = [allele1, allele2, phased]
        a1, a2 = gt[0], gt[1]

        if a1 == -1 or a2 == -1:
            gt_class = "MISSING"
        elif a1 == 0 and a2 == 0:
            gt_class = "HOM_REF"
        elif a1 != a2:
            gt_class = "HET"
        elif a1 > 0 and a2 > 0:
            gt_class = "HOM_ALT"
        else:
            gt_class = "OTHER"

        rows.append([
            var.CHROM,
            var.POS,
            var.ID,
            var.REF,
            ",".join(var.ALT),
            str(sample).strip(),
            gt_class
        ])

vcf.close()

df = pd.DataFrame(
    rows,
    columns=["CHROM", "POS", "ID", "REF", "ALT", "SAMPLE", "GENOTYPE_CLASS"]
)

pivot = df.pivot(
    index=["CHROM", "POS", "ID", "REF", "ALT"],
    columns="SAMPLE",
    values="GENOTYPE_CLASS"
)

pivot = pivot.sort_index(axis=1)
pivot = pivot.fillna("MISSING")

pivot.to_csv("per_variant_samples_genotype_class.tsv", sep="\t")

gi_info = pd.read_excel("pop_info.xlsx")
gi_info = gi_info.drop(columns=["IID", "Order", "Demography"])

counts = {
    str(s).strip(): defaultdict(int)
    for s in samples
}

for var in vcf:
    for sample, gt in zip(samples, var.genotypes):
        sample = str(sample).strip()
        a1, a2 = gt[0], gt[1]

        if a1 == -1 or a2 == -1:
            counts[sample]["MISSING"] += 1
        elif a1 == 0 and a2 == 0:
            counts[sample]["HOM_REF"] += 1
        elif a1 != a2:
            counts[sample]["HET"] += 1
        elif a1 > 0 and a2 > 0:
            counts[sample]["HOM_ALT"] += 1

vcf.close()

geno_counts = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
geno_counts = geno_counts.reset_index().rename(columns={"index": "SAMPLE"})

for col in ["HOM_REF", "HET", "HOM_ALT", "MISSING"]:
    if col not in geno_counts.columns:
        geno_counts[col] = 0

geno_counts = geno_counts[["SAMPLE", "HOM_REF", "HET", "HOM_ALT", "MISSING"]]

merged = geno_counts.merge(
    gi_info,
    left_on="SAMPLE",
    right_on="FID",
    how="left"
).drop(columns=["FID"])

merged.to_csv(
    "per_sample_clinvar_p_lp_genotype_counts_annotated.tsv",
    sep="\t",
    index=False
)