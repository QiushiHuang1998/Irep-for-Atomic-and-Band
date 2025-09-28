from pymatgen.core.structure import Structure, Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, PointGroupAnalyzer
from collections import OrderedDict
import pandas as pd
import sys
import re
import os
# read POSCAR
original_stdout = sys.stdout

# 指定输出文件的路径
output_file_path = "output.txt"

try:

# 打开输出文件，将标准输出流重定向到文件
    with open(output_file_path, "w") as file:
        sys.stdout = file  
        current_directory = os.getcwd()
        poscar_file = os.path.join(current_directory, "POSCAR")

        with open(poscar_file, "r") as file:
            lines = file.readlines()
        elements = lines[5].split()
        atom_counts = list(map(int, lines[6].split()))
        element_list = []
        for i in range(len(elements)):
            element_list.extend([elements[i]] * atom_counts[i])

        structure = Structure.from_file(poscar_file)
        filename = "./data"


        def convert_letter_to_number(letter):
            if len(letter) != 1 or not letter.isalpha() or not letter.islower():
                return None
            else:
                return ord(letter) - ord('a')
        # spd to sitegroup irrep
        def search_character(filename, character):
            with open(filename, 'r') as file:
                lines = file.readlines()

                for i, line in enumerate(lines):
                    if line.strip() == character:
                        start_index = i
                        break

                if start_index is not None:
                    end_index = min(start_index + 6, len(lines))
                    result = lines[start_index:end_index]
                    return result
                else:
                    return None
        # sitegroup irrep to spacegroup irrep
        def search_character_in_dataframe(df, row_char, col_char):
            row_indices = df[df.astype(str).apply(lambda x: row_char in ' '.join(x), axis=1)].index
            col_indices = df.columns[df.columns.astype(str).str.contains(col_char) & ~df.columns.astype(str).str.contains("Band")]
        # print("dg",row_indices,col_indices)   
            if len(row_indices) > 0 and len(col_indices) > 0:
                result = df.loc[row_indices[0], col_indices[0]]
                return result

            return None


        # 分析对称性和空间群
        analyzer = SpacegroupAnalyzer(structure)
        pointgroup = analyzer.get_point_group_symbol()
        print("Point Group: ", pointgroup)
        symstructure=analyzer.get_symmetrized_structure()
        print(symstructure)
        spacegroup = analyzer.get_space_group_symbol()
        wyckoff_sites = analyzer.get_symmetrized_structure().equivalent_sites
        site_nums = analyzer.get_symmetrized_structure().wyckoff_symbols
        #print(site_nums)
        site_groups = analyzer.get_symmetry_dataset()["site_symmetry_symbols"]
        #print(site_groups)
        spacegroup_number = analyzer.get_space_group_number()
        site_numses =analyzer.get_symmetry_dataset()["wyckoffs"]
        comb = [f"{site_groups[i]} {site_numses[i]}" for i in range(len(element_list))]
        #print(comb)
        unique_elements = list(OrderedDict.fromkeys(comb))
        #print(unique_elements)
        # 获取每个Site的点群符号
        #print("空间群编号:", spacegroup_number)
        #print("空间群: ", spacegroup)
        l=0
        m=0
        for sites in wyckoff_sites:
            l +=1
            m += len(sites)
            print("Site Group: ", site_nums[l-1],element_list[m-1])
            for site in sites:
                print("Site: ", [site.a,site.b,site.c])
        
        a = spacegroup_number
        df0 = pd.read_csv(f"./tabletest/{a}-a.txt", delimiter='\t')
        print("\n",spacegroup_number,"Spacegroup-High Symmetry Kpoints：",df0.iloc[2:, 0].values)


        row_char = "0,0,0"  #Kpoints Input
        z=0
        for unique_element in unique_elements:
            z += 1
            parts = unique_element.split()
            letter = str(parts[1])
            letters = re.sub(r'[^a-zA-Z]', '', letter)
            number = convert_letter_to_number(letters)
            
            b = letters
            df = pd.read_csv(f"/public23/home/sca1074/code/xuxun/xuxun2/tabletest/{a}-{b}.txt", delimiter='\t')
            new_string = re.sub(r'\.', '', parts[0])
            #print(new_string)
            character = str(new_string)
            result = search_character(filename, character)
            if result is not None:
                print("\n",letter,"site Atomic orbitals occupy the irreducible representations of the site group：")
                formatted_result = '\n'.join(['\t'.join(line.strip().split()) for line in result])
                print(formatted_result) 
                lines = formatted_result.strip().split('\n')
                print("\n",letter,"site Atomic orbitals occupy the irreducible representations of the energy band.：")
                for line in lines[2:]:
                    columns = line.strip().split('\t')
                    orbital = str(columns[0])
                    result_list = []
                    for column in columns[1:]:
                        col_char = str(column)
                    # print("debug",row_char,col_char)
                        result = search_character_in_dataframe(df, row_char, col_char)
                        result_list.append(result)
                    print(orbital, "\t".join(result_list))
            else:
                print("ero1")
except Exception as e:
    print(current_directory)

finally:
    # 恢复标准输出流
    sys.stdout = original_stdout
    print(current_directory)
# 输出已保存到文件
print(f"Output has been saved to {output_file_path}")
