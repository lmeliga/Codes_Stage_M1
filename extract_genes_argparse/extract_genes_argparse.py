from Bio import SeqIO #Ce module SeqIO de Biopython permet de lire et d'écrire des fichiers de séquence dans différents formats
from BCBio import GFF # Ce module "GFF" sert à lire et manipuler les fichiers d'annotation .gff 
from Bio.Seq import Seq
import argparse 

def main():
    parser = argparse.ArgumentParser(description='Use fasta file and gff file to "extract genes" to find which part of the sequence match to which gene')
    parser.add_argument('fasta_file', type=str, help = 'fasta file which contains a sequence')
    parser.add_argument('gff_annotations', type=str, help = 'GFF file with the complete annotations')
    args = parser.parse_args() # analyse + stockage des valeurs entrées dans args 
    fasta_file = args.fasta_file # on assigne les valeurs des arguments à fasta_file
    gff_annotations = args.gff_annotations # on assigne les valeurs des arguments à gff_annotations

    seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")) # lecture de la séquence du fichier fasta et stockage dans un dictionnaire (ici, "seq_records") où la clé est l'identifiant de séquence et la valeur est l'objet seq_record correspondant 

    with open(gff_annotations) as gff: # ouvrir le fichier .gff et le lire 
        for rec in GFF.parse(gff, base_dict=seq_records): #parse les annotations et les associe aux séquences correspondantes dans seq_records. "base_dict=seq_records" indique au parseur d'utiliser les séquences du dictionnaire seq_records comme base.
            seq_record = rec # stocke chaque enregistrement annoté dans seq_record

    for feature in seq_record.features: # itération sur toutes les caractéristiques annotées dans la séquence. L'attribut "features" provient de SeqRecord
        # type, location.start, location.end et strand sont des attributs de SeqFeature
        if feature.type == "gene": # On ne s'intéresse qu'aux caractéristiques de type "gene" 
            start = feature.location.start # position du début du gène 
            end = feature.location.end # position de la fin du gène
            strand = feature.strand # indique le type de brin
            # gene_seq = seq_record.seq[start:end] # extraction de la sequence entre le début et la fin 
            gene_seq = seq_record.seq[start:end]

            if strand == -1:
                gene_seq = Seq.reverse_complement(gene_seq) #si le gène est sur le brin antisens (strand == -1), on prend la séquence complémentaire 

            gene_id = feature.qualifiers.get('ID', ['unknown'])[0] #récupère l'ID du gène à partir du fichier gff d'annotations (unknown si pas d'ID)
            print(f"Gene {gene_id}: {gene_seq}\n") #affiche l'ID gu gène + séquence 

if __name__ == "__main__": # Vérification de l'exécution du script 
    main() # Appel de la fonction principale main()
