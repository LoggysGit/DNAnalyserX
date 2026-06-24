def format_mutation_results(raw_res):
    if not raw_res: return []

    res = []
    # Put first element
    temp_mut = [raw_res[0][0], raw_res[0][1], raw_res[0][2], raw_res[0][3]]

    for i in range(1, len(raw_res)):
        mutation = raw_res[i]

        if mutation[0] == temp_mut[0] + len(temp_mut[2]) and mutation[1] == temp_mut[1]:
            # Don't copy gaps (dots)
            if mutation[2] != ".": temp_mut[2] += mutation[2]
            if mutation[3] != ".": temp_mut[3] += mutation[3]
        else:
            # Single -> Multiple
            if temp_mut[1] == "SNP" and (len(temp_mut[2]) > 1 or len(temp_mut[3]) > 1): temp_mut[1] = "MNP"

            res.append(temp_mut)

            temp_mut = [mutation[0], mutation[1], mutation[2], mutation[3]]

    if temp_mut[1] == "SNP" and (len(temp_mut[2]) > 1 or len(temp_mut[3]) > 1): temp_mut[1] = "MNP"
    res.append(temp_mut)

    return res

raw_startmerged = [
     [65100, "SNP", "A", "C"],
     [65101, "SNP", "A", "T"],
     [65102, "SNP", "G", "A"],
     [66300, "SNP", "A", "C"],
     [68320, "Deletion", "A", "."],
     [70100, "Insertion", ".", "C"],
     [71320, "Deletion", "T", "."],
     [71321, "Deletion", "G", "."],
     [72100, "SNP", "A", "C"],
]

raw_lastmerged = [
     [65100, "SNP", "A", "C"],
     [66300, "SNP", "T", "C"],
     [68320, "Deletion", "A", "."],
     [70100, "Insertion", ".", "C"],
     [72100, "SNP", "A", "C"],
     [72101, "SNP", "G", "T"],
]

raw_lasttwomerged = [
     [65100, "SNP", "A", "C"],
     [66300, "SNP", "T", "C"],
     [68320, "Deletion", "A", "."],
     [70100, "Insertion", ".", "C"],
     [70101, "Insertion", ".", "A"],
     [72100, "SNP", "A", "C"],
     [72101, "SNP", "G", "T"],
]

print(format_mutation_results(raw_startmerged))
print(format_mutation_results(raw_lastmerged))
print(format_mutation_results(raw_lasttwomerged))