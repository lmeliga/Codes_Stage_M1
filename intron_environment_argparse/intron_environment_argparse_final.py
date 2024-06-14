import re # Importation du module pour gérer les expressions régulières
import os # importation du module os pour les opérations système 
import argparse # pour la gestion des arguments de la ligne de commande 

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

def gff2_parser(gff_file2, start_position, end_position, distance):
    upstream_genes = []
    downstream_genes = []

    upstream_range = (start_position - distance, start_position)
    downstream_range = (end_position, end_position + distance)

    with open(gff_file2, 'r') as gff:
        for line in gff:
            colonne = re.split(r'[ =\t;]+', line)
            if len(colonne) < 5 or colonne[2] != 'CDS':
                continue
            gene_start = int(colonne[3])
            gene_end = int(colonne[4])

            if upstream_range[0] <= gene_start <= upstream_range[1] or upstream_range[0] <= gene_end <= upstream_range[1]:
                upstream_genes.append(line.strip())
            if downstream_range[0] <= gene_start <= downstream_range[1] or downstream_range[0] <= gene_end <= downstream_range[1]:
                downstream_genes.append(line.strip())

        return upstream_genes, downstream_genes

def main(): # boucle principale 
    parser = argparse.ArgumentParser(description='Use GFF files to find upstream and downstream genes.') # Création d'un objet "argumentparser" pour gérer les arguments de la ligne de commande  
    parser.add_argument('gff_file1', type=str, help='GFF file with the data of the insertion')# ajout d'un argument pour le premier gff ( = celui avec les informations sur l'insertion)
    parser.add_argument('gff_file2', type=str, help='GFF file with the complete annotations') # ajout d'un argument pour le deuxième gff ( = celui avec les informations sur l'annotation complète correspondante au gff1)
    parser.add_argument('distance', type =int, help = 'Maximum distance upstream and downstream the intron at which we will look for genes')

    args = parser.parse_args() # analyse + stockage des valeurs entrées dans args 

    gff_file1 = args.gff_file1 # on assigne les valeurs des arguments à gff_file1
    gff_file2 = args.gff_file2 # idem pour gff_file2
    distance = args.distance # distance maximale à laquelle on va chercher des gènes en amont et en aval de l'intron
    upstream_counter = 0
    downstream_counter = 0

################################### VERIFICATION/CONTROLE ##########################################

    if not (gff_file1.endswith(".gff") and os.path.isfile(gff_file1)):
        print("Le fichier gff_file1 n'est pas valide ou ne porte pas l'extension .gff")
        return

    if not (gff_file2.endswith(".gff") and os.path.isfile(gff_file2)):
        print("Le fichier gff_file2 n'est pas valide ou ne porte pas l'extension .gff")
        return

    with open(gff_file1, 'r') as gff:
        column3 = [line.split("\t")[2] for line in gff if not line.startswith("#") and len(line.split("\t")) > 2]

    if "insertion" not in column3:
        print("Aucune entrée 'insertion' n'a été trouvée dans le gff_file1 fourni.")
        return

    with open(gff_file2, 'r') as gff:
        if not any(line.startswith("#") for line in gff):
            print("Le fichier gff_file2 ne contient pas d'annotations.")
            return
################################### VERIFICATION/CONTROLE ##########################################

    start_position, end_position = gff1_parser(gff_file1)

    if start_position is not None and end_position is not None:
        print("Position de début pour waldern_only21nt5 : ", start_position)
        print("Position de fin pour waldern_3end : ", end_position)
        print()
        print()
        print()
        

        genes_upstream, genes_downstream = gff2_parser(gff_file2, start_position, end_position, distance)

        if genes_upstream is not None and genes_downstream is not None:
            print(f"Gènes en amont (dans la plage start_position - {distance} à start_position) :")
            for gene in genes_upstream:
                upstream_counter += 1
                print(gene)
                

            # Les print suivants servent à marquer séparation entre les gènes en amont et en aval
            print()
            print()
            print()

            print(f"Gènes en aval (dans la plage end_position à end_position + {distance}) :")
            for gene in genes_downstream:
                downstream_counter += 1
                print(gene)
            print()
            print()
            print()
            
            print(f"Il y a {upstream_counter} gènes en AMONT.")
            print(f"Il y a {downstream_counter} gènes en AVAL.")
        else:
            print("Les gènes n'ont pas pu être extraits correctement.")
    else:
        print("Les positions n'ont pas pu être extraites correctement.")

if __name__ == "__main__": # Vérification de l'exécution du script 
    main() # Appel de la fonction principale main()
