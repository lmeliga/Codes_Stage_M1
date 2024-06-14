import re
import os 
import argparse

def gff1_parser(gff_file1):
    start_position = None
    end_position = None
    with open(gff_file1, 'r') as gff:
        for line in gff:
            colonne = line.split()
            if len(colonne) < 5:
                 continue
            if "waldern_only21nt5" in colonne[2]:
                        start_position = min(int(colonne[3]), int(colonne[4]))
            elif "waldern_3end" in colonne[2]:
                        end_position = max(int(colonne[3]), int(colonne[4]))
        if start_position is not None and end_position is not None:
            return start_position, end_position
        else:
            print("Les lignes waldern_only21nt5 ou waldern_3end n'ont pas été trouvées correctement.")
            return None, None

def gff2_parser(gff_file2, start_position, end_position):
    btw_genes = [] # On crée la liste des genes qui sont dans l'intron. Au départ, elle est vide

    if start_position > end_position :
        range_limit = (end_position, start_position)
    else : 
        range_limit = (start_position, end_position)

    with open(gff_file2, 'r') as gff:
        for line in gff:
            colonne = re.split(r'[ =\t;]+', line)
            if len(colonne) < 5 or colonne[2] != 'gene':
                    continue  # Ignore la ligne si jamais il n'y a pas "gene" 
            gene_start = int(colonne[3])
            gene_end = int(colonne[4])

            if range_limit[0] <= gene_start <= range_limit[1] and range_limit[0] <= gene_end <= range_limit[1]:
                    btw_genes.append(line.strip()) # si gene_start et gene_end sont entre les deux bornes, dans ce cas on ajoute la ligne entière à la liste

        return btw_genes

def gene_finder_ID(gff_file, gene_wanted): # la fonction prend en arguments le fichier gff avec les annotations (le GCA) et le gène recherché (attention à l'orthographe)
    count_genes = 0 # création et initialisation d'un compteur pour compter le nombre de gènes recherchés présents dans le gff au total
    gene = [] # création d'une liste vide qui contiendra les lignes des gènes visés
    with open(gff_file, 'r') as gff : # on ouvre le fichier gff en question (on le renomme "gff" juste dans la fonction pour faciliter l'écriture)
         for line in gff: # Pour chaque ligne du fichier, on effectue la commande en dessous
            if not line.startswith("#") : # Si la ligne ne commence PAS par un #...
                colonne = line.strip().split("\t") # On définit une variable colonne. On divise la ligne en colonnes en utilisant la tabulation en guise de séparateur et on retire tout potentiel espace en début et en fin de ligne 
                if len(colonne) < 9 or colonne[2] != 'CDS': # Si la ligne possède moins de 9 colonnes ou qu'il n'y a pas écrit "CDS" à la troisième colonne (normalement il y a écrit "gene" ou "CDS")
                    continue # Si une des deux conditions précédentes n'est pas respectée, on passe à la ligne suivante
                qualifiers = "\t".join(colonne[8:]) # On fusionne toutes les colonnes à partir de la 9ème jusqu'à la dernière (peu importe le nombre) en une seule chaine de caractères en utilisant la tabulation comme séparateur 
                second_part = re.split(r'[ =\t:;()]+', qualifiers) # définit la "seconde partie" de la ligne, et divise uniquement cette dernière (grâce à qualifiers) selon les séparateurs qui sont l'espace, =, la tabulation, :, ;, (), et les nombres de 0 à 9. Le "+" indique qu'il peut y avoir plusieurs de chacun de ces séparateurs
                if gene_wanted in second_part: # Si le gène recherché est dans la "seconde partie"...
                    gene.append(line.strip()) # ...on ajoute la ligne dans la liste nommée "gene"
                    count_genes += 1 # ... et on incrémente le compteur de 1
    return gene, count_genes


def main():
    
    parser = argparse.ArgumentParser(description='Use gff files to find all of the genes inside the intron itself')
    parser.add_argument('gff_file1', type=str, help='Gff file with informations about insertion')
    parser.add_argument('gff_file2', type=str, help='Gff file with annotations')
    args = parser.parse_args()
    gff_file1 = args.gff_file1
    gff_file2 = args.gff_file2

    ################################### VERIFICATION/CONTROLE ##########################################

    if not (gff_file1.endswith(".gff") and os.path.isfile(gff_file1)):
        print("Le premier fichier gff n'est pas valide ou ne porte pas l'extension .gff")
        return

    if not (gff_file2.endswith(".gff") and os.path.isfile(gff_file2)):
        print("Le deuxième fichier gff n'est pas valide ou ne porte pas l'extension .gff")
        return

    with open(gff_file1, 'r') as gff:
        column3 = [line.split("\t")[2] for line in gff if not line.startswith("#") and len(line.split("\t")) > 2]

    if "insertion" not in column3:
        print("Aucune entrée 'insertion' n'a été trouvée dans le gff_file1 fourni.")
        return

    with open(gff_file2, 'r') as gff:
        if not any(line.startswith("#") for line in gff):
            print("Le deuxième gff ne contient pas d'annotations.")
            return

    ################################### VERIFICATION/CONTROLE ##########################################

    start_position, end_position = gff1_parser(gff_file1)

    genes_id = []

    if start_position is not None and end_position is not None:
        print("L'intron commence en position : ", start_position)
        print("L'intron se termine en position : ", end_position)
        print("Détermination des gènes entre ces deux bornes...")

        genes_in_between = gff2_parser(gff_file2, start_position, end_position)

        if genes_in_between is not None :
            print(" Les gènes présents dans l'intron sont : ")
            for gene in genes_in_between:
                colonne = gene.split()
                gene_qualifiers = colonne[8]
                try:
                    gene_id = re.search('ID=gene-(.+?);', gene_qualifiers).group(1)
                except AttributeError:
                    gene_id = 'ID du gène non trouvé'
                print(gene_id)
                genes_id.append(gene_id)
            
            for gen_id in genes_id:
                gen,number = gene_finder_ID(gff_file2,gene_id)
                upgraded_genes = "\n\n".join(gen)
                print("Gènes trouvés :")
                print(upgraded_genes)
                print("\nNombre de gènes trouvés :", number)
        else:
            print("Les gènes n'ont pas pu être extraits correctement.")
    else:
        print("Les positions n'ont pas pu être extraites correctement.")

if __name__ == "__main__":
    main()