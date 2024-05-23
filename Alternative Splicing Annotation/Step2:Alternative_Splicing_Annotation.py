def process_file(filename, output_filename):
    with open(filename, 'r', encoding='utf-8') as file, open(output_filename, 'w', encoding='utf-8') as outfile:
        lines = file.readlines()
        outfile.write("chrom" + "\t" + "Cluster" + "\t" + "ASevent" + "\n")
        for i in range(3, len(lines), 5):  
            parts = lines[i].strip().split(':')
            if len(parts) < 3:
                continue  
            key = parts[0].strip()     
            key = key.replace("(", "").replace(")", "").replace("'", "")
            key_parts = key.split(',')
            chrom = key_parts[0].strip()
            cluster = key_parts[1].strip()
            anno = parts[2].strip()
            outfile.write(f"{chrom}\t{cluster}\t{anno}\n")




if __name__ == "__main__":
    input_filename = 'annotated_position.txt' 
    output_filename = 'annotated_data.txt' 

    process_file(input_filename, output_filename)
    print(f"Processed data has been saved to {output_filename}")





