from Bio import SeqIO # Lire les séquences fasta
from Bio.Seq import Seq
import re 
import argparse

def main() : 
    parser = argparse.ArgumentParser(description='Use fasta file and gff file to find the original STOP codon that existed in the sequence before the insertion of the intron')
    parser.add_argument('fasta_file', type=str, help = 'fasta file which contains a sequence')
    parser.add_argument('gff_insertion', type=str, help = 'GFF file with the informations on the insertion in the intron')
    parser.add_argument('motif', type=str, help='Pattern that you are looking for in the sequence')
    args = parser.parse_args() # analyse + stockage des valeurs entrées dans args 
    fasta_file = args.fasta_file # on assigne les valeurs des arguments à fasta_file
    gff_insertion = args.gff_insertion # on assigne les valeurs des arguments à gff_insertion
    motif = args.motif # on assigne les valeurs des arguments à motif 


    # Création et initialisation de l'ensemble des variables qui serviront durant le programme 
    wonly21_debut = None 
    w3end_debut = None 
    start_only21 = None 
    end_only21 = None
    start_3end = None 
    sense_strand = True # Ce booléen servira à diviser le programme en deux. Voir plus bas. Par défaut, on considère que la séquence est lue dans le sens 5' -> 3' et donc se trouve sur le brin sens.
    seq_records = list(SeqIO.parse(fasta_file, "fasta")) # Chaque séquence du fasta est enregistrée dans seq_records. Dans le contexte actuel, il n'y en a qu'une 


    l_STOP_before = [] # Création et initialisation de la liste qui contiendra les nucléotides du codon STOP initial  

    with open(gff_insertion, 'r') as gff : # Le gff est ouvert une première fois pour extraire les positions du début de only21 et de 3end 
        for line in gff:
            colonne = line.split("\t")
            if colonne[2] == "waldern_only21nt5":
                wonly21_debut = min(int(colonne[3]), int(colonne[4]))
            if colonne[2] == "waldern_3end":
                w3end_debut = min(int(colonne[3]), int(colonne[4]))

    if wonly21_debut > w3end_debut : # Les deux positions sont comparées pour savoir dans quel brin on se trouve, et le booléen est modifié en conséquence 
        sense_strand = False 
    else : 
        sense_strand = True

    # A partir d'ici, le programme est divisé en deux et fonctionnera différemment en fonction de si on se trouve sur le brin sens (sense_strand = True) ou antisens (sense_strand = False).

    if sense_strand : # Nous sommes dans le cas où nous sommes dans le brin sens. Cela est précisé à l'utilisateur via la ligne suivante 
        print("La séquence se trouve sur le brin sens.")

        # Ces variables ne serviront QUE dans cette situation
        start_only21 = None 
        end_final = None 

        with open(gff_insertion, 'r') as gff : # Le gff est à nouveau ouvert pour extraire les positions du début du only21 et de la fin du 3end, qui correspond exactement à la fin de l'intron. 
            for line in gff:
                colonne = line.split("\t")
                if colonne[2] == "waldern_only21nt5":
                    start_only21 = min(int(colonne[3]), int(colonne[4]))
                elif colonne[2] == "waldern_3end":
                    end_final = max(int(colonne[3]), int(colonne[4]))
                if start_only21 is not None and end_final is not None : 
                    break # Quand les deux positions ont été trouvées, on sort de la boucle 

        if start_only21 is None or end_final is None : 
            raise ValueError("Les annotations waldern_only21nt5 ou waldern_3end sont manquantes dans le fichier GFF que vous avez entré.") # Si au moins une des deux positions n'a pas été trouvée, l'erreur s'affiche 


        for seq_record in seq_records: # On s'intéresse à chaque séquence enregistrée dans seq_records. Ici, il n'y en a qu'une donc la boucle n'itérera qu'une fois.
            sequence = str(seq_record.seq) # La séquence est transformée en chaine de caractères
            print(f"Recherche du motif dans la sequence : {seq_record.id}")

            for match in re.finditer(motif, sequence): # Utilise "re.finditer" pour trouver toutes les occurences du motif dans la séquence 
                start = match.start() # Obtient la position du début du motif trouvé
                end = match.end() # Obtient la position de la fin du motif trouvé 
                if end >= start_only21 and end <= end_final: # Normalement, la position de la fin du motif d'intérêt (donc du début de l'intron) entre le début de only21 et la fin de 3end 
                    print(f"Le motif représentant le début de l'intron se situe entre {start} et {start + len(motif)} : {sequence[start:start+len(motif)]}")
                    if start >= 2 : # Simple vérification pour nous assurer qu'il y a assez de nucléotides avant le start 
                        # On récupère les deux nucléotides avant le début du motif 
                        nucleotide_minus_two = sequence[start - 2]
                        nucleotide_minus_one = sequence[start - 1]
                    else : 
                        raise ValueError("Impossible d'extraire les nucléotides avant le motif : indices négatifs")
                    
                    if end_final < len(sequence): #éviter indice hors limite 
                        nucleotide_plus_one = sequence[end_final] 
                    else : 
                        raise ValueError("indice hors limite")
                    
                    l_STOP_before.extend([nucleotide_minus_two, nucleotide_minus_one, nucleotide_plus_one]) # Ajout des nucléotides trouvés dans la liste de sorte à former le codon
                    break # Le motif d'intérêt a été trouvé + les nucléotides ont été extraits + le codon est formé => on sort de la boucle.

        print("Nucléotides avant l'insertion de l'intron:", l_STOP_before) # On affiche le codon

    else : # Nous sommes dans le cas où la séquence fasta se trouve sur le brin antisens 
        print("La séquence se trouve sur le brin antisens.") 
        new_motif = str(Seq(motif).reverse_complement()) # Convertir en objet Seq, puis utiliser reverse_complement, et convertir en string
        
        with open(gff_insertion, 'r') as gff : 
            for line in gff:
                colonne = line.split("\t")
                if colonne[2] == "waldern_only21nt5":
                    start_only21 = min(int(colonne[3]), int(colonne[4]))
                    end_only21 = max(int(colonne[3]), int(colonne[4]))
                elif colonne[2] == "waldern_3end":
                    start_3end = min(int(colonne[3]), int(colonne[4]))
                if end_only21 is not None and start_3end is not None and start_only21 is not None : 
                    break 

        if end_only21 is None or start_3end is None or start_only21 is None : 
            raise ValueError("Les annotations waldern_only21nt5 ou waldern_3end sont manquantes dans le fichier GFF que vous avez entré.") 

        for seq_record in seq_records:
            sequence = str(seq_record.seq)
            print(f"Recherche du motif dans la sequence : {seq_record.id}")

            for match in re.finditer(new_motif, sequence):
                start = match.start()
                end = match.end()
                if end >= start_only21 and end <= end_only21: 
                    print(f"Le motif représentant le début de l'intron se situe entre {start} et {start + len(new_motif)} : {sequence[start:start+len(new_motif)]}")
                    # En terme de raisonnement, il faut raisonner à l'envers, puisque le brin est lu à l'envers du sens habituel.
                    if end >= 2 : 
                        nucleotide_plus_two = sequence[end + 3]
                        nucleotide_plus_one = sequence[end + 1]
                    else : 
                        raise ValueError("Impossible d'extraire les nucléotides avant le motif : indices négatifs")
                    
                    if start_3end < len(sequence): #éviter indice hors limite 
                        nucleotide_minus_one = sequence[start_3end - 1] 
                    else : 
                        raise ValueError("indice hors limite")
                    
                    l_STOP_before.extend([nucleotide_minus_one, nucleotide_plus_one, nucleotide_plus_two])
                    break

        print("Nucléotides avant l'insertion de l'intron:", l_STOP_before)

if __name__ == "__main__": # Vérification de l'exécution du script 
    main() # Appel de la fonction principale main()
