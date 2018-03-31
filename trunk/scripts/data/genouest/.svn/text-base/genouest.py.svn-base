import pickle

def build_geneoust():

    def parse_geneoust(data_type):
        db = {}

        bacteria_source_filename_split = open("bacteria_" + data_type + ".txt").read().split('>')
        archaea_source_filename_split = open("archaea_" + data_type + ".txt").read().split('>')
        
        source_filename_split = bacteria_source_filename_split + archaea_source_filename_split
            
        for a in source_filename_split:
            b = a.split('|')
            if len(b) == 3:
                tax_id = b[1]
                tax_name = ' '.join(b[2].split('\n')[0].split(' ')[:-1])[2:]
                coords = b[2].split('\n')[0].split(' ')[-1]
                coords = coords.split('-')
                coords = [int(coords[0]), int(coords[1])]
                
                seq = b[2].split('\n')[1]
                if tax_id not in db.keys():
                    db[tax_id] = [[], []]
                if tax_name not in db[tax_id][0]:
                    db[tax_id][0].append(tax_name)
                db[tax_id][1].append([coords, seq])
        wfile = open("geneoust_" + data_type + "_db.txt", "w")
        pickle.dump(db, wfile)
        wfile.close()

    parse_geneoust("repeats")
    parse_geneoust("spacers")

def build_array_db():    
    repeats_db = pickle.load(open("geneoust_repeats_db.txt"))
    spacers_db = pickle.load(open("geneoust_spacers_db.txt"))

    repeats_db_keys = repeats_db.keys()
    spacers_db_keys = spacers_db.keys()

    array_dict = {}

    for a in repeats_db_keys:
        if a in spacers_db_keys:
            array_dict[a] = [[repeats_db[a][0], spacers_db[a][1]], []]

            for b in repeats_db[a][1]:
                array_dict[a][1].append([b[0][0], [b, 'repeat']])
            for b in spacers_db[a][1]:
                array_dict[a][1].append([b[0][0], [b, 'spacer']])

            array_dict[a][1].sort()
            
    wfile = open("geneoust_array_db.txt", "w")
    pickle.dump(array_dict, wfile)
    wfile.close()
