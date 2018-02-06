import re

class glymmer_module:

    def __init__(self, parameters, max_mass):
        self.parameters = parameters
        self.max_mass = max_mass
        self.regex_rules = []
        self.regex_rules.append('\(Ma3\(Ma2Ma3\(Ma6\)Ma6') #rule 2 DONE
        self.regex_rules.append('\(Ma2Ma') #rule 1 DONE
        self.regex_rules.append('[ab][0-9]\(GNb[0-9]\)\(.*\)') #rule 5, bisecting glcnac DONE
        self.regex_rules.append('\(Ma3\(Ma3\(Ma6\)Ma6\)Mb4') #rule 8 DONE
        self.regex_rules.append('\(Ab4GN') #rule 15 DONE
        self.regex_rules.append('\(.*Ab4GN') #rule 24 DONE
        self.regex_rules.append('\(Ab4GN') #rule 19 DONE
        self.regex_rules.append('\(Ab4GN') #rule 17 DONE
        self.regex_rules.append('\(Ab4GNb3Ab') #rule 16 DONE
        self.regex_rules.append('[ab][0-9]\(GNb[0-9]\)\(.*\)') #rule 18
        self.regex_rules.append('\(GNb2[\(.*\)|\)]*Ma3\(Ma6\)Mb4') #rule 9 DONE
        self.regex_rules.append('\(GNb2Ma3') # rule 11 DONE
        self.regex_rules.append('\(GNb2Ma6') # rule 12 DONE
        self.regex_rules.append('\(Ab4GN') # rule 13 NOT DONE (CONSTRAINT)
        self.regex_rules.append('\(Ma3\(Ma6\)Ma6') # rule 5 DONE
        self.regex_rules.append('\(Ma6Ma6') # rule 6 DONE
        self.regex_rules.append('\(GN') # rule 14 and rule 18 NOT DONE (CONSTRAINT)
        self.regex_rules.append('GNb4GN') # rule 7 DONE
        self.regex_rules.append('\(Fa2Ab') # rule 25 and rule 26
        self.regex_rules.append('\(Ab4GN') #rule 19
        self.regex_rules.append('Ab3GNb') # rule 20
        self.regex_rules.append('\(.*Ab4GNb') # rule 21 and 24
        self.regex_rules.append('\(Ab3GNb') # rule 22
        self.regex_rules.append('\(Ab4GNb') # rule 23
        self.regex_rules.append('\(NNa3Ab3AN')# rule 27
        self.regex_rules.append('GNb2[\(.*\)|\)]*Ma3') # rule 10


    def generate_network(self, substrates):
        count = 0
        f = open("reactions.txt", "w")
        seen_substrates = set()
        for substrate in substrates:
            if substrate not in seen_substrates:
                seen_substrates.add(substrate)
                reactions = []
                print "loop iteration", count
                new_substrates = []
                substrate = "(" + substrate
                if "Ga3" not in substrate:
                    search_string = ""
                    #Man1 - specific structure
                    #one = True
                    if re.findall(self.regex_rules[0], substrate) != []:
                        ind = re.search(self.regex_rules[0], substrate).span()
                        search_string = substrate[ind[0]:ind[1]]
                        print "rule 2"
                        #one = False
                    #Man1 - general removal of Ma2
                    else:
                        if re.findall(self.regex_rules[1], substrate) != []:
                            ind = re.search(self.regex_rules[1], substrate).span()
                            search_string = substrate[ind[0]:ind[1]]
                            print "rule 1"
                    # actually removes Ma2
                    if search_string != "":
                        starts_one = [m.start() for m in re.finditer(self.regex_rules[1], substrate)]
                        for start in starts_one:
                            prod = substrate[:start+1] + substrate[start+4:]
                            reactions.append("rule 1 " + substrate + " " + prod[1:])
                            new_substrates.append(prod[1:])
                        starts_two = [m.start() for m in re.finditer(self.regex_rules[0], substrate)]
                        for start in starts_two:
                            prod = substrate[:start+6] + substrate[start+9:]
                            reactions.append("rule 2 " + substrate + " " + prod[1:])
                            new_substrates.append(prod[1:])
                #GnT1 action
                if re.findall(self.regex_rules[3], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[3], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = "(GNb2" + substrate[start_split:]
                        reactions.append("rule 8 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                # GnT11 action
                # check for no bisecting glcnac
                if re.findall(self.regex_rules[2], substrate) == []:
                    if re.findall(self.regex_rules[10], substrate) != []:
                        ends = [m.end() for m in re.finditer(self.regex_rules[10], substrate)]
                        for end in ends:
                            end_split = end - 7
                            prod = substrate[:end_split] + "GNb2" + substrate[end_split:]
                            reactions.append("rule 9 " + substrate + " " + prod[1:])
                            new_substrates.append(prod[1:])
                    if re.findall(self.regex_rules[11], substrate) != []:
                        ends = [m.end() for m in re.finditer(self.regex_rules[11], substrate)]
                        for end in ends:
                            end_split = end - 3
                            prod = substrate[:end_split] + "(GNb4)" + substrate[end_split:]
                            reactions.append("rule 11 " + substrate + " " + prod[1:])
                            new_substrates.append(prod[1:])
                    if re.findall(self.regex_rules[12], substrate) != []:
                        ends = [m.end() for m in re.finditer(self.regex_rules[12], substrate)]
                        for end in ends:
                            end_split = end - 3
                            prod = substrate[:end_split] + "(GNb6)" + substrate[end_split:]
                            reactions.append("rule 12 " + substrate + " " + prod[1:])
                            new_substrates.append(prod[1:])
                    if re.findall('\(GNb2[\(.*\)|\)]*Ma3', substrate) != []:
                        if re.findall(self.regex_rules[14], substrate) != []:
                            inds = [(m.start(),m.end()) for m in re.finditer(self.regex_rules[14], substrate)]
                            for ind in inds:
                                prod = substrate[:ind[0]] + "(Ma6Ma6" + substrate[ind[1]:]
                                reactions.append("rule 5 " + substrate + " " + prod[1:])
                                new_substrates.append(prod[1:])
                        if re.findall(self.regex_rules[15], substrate) != []:
                            #rule 6
                            inds = [(m.start(),m.end()) for m in re.finditer(self.regex_rules[15], substrate)]
                            for ind in inds:
                                prod = substrate[:ind[0]] + "(Ma6" + substrate[ind[1]:]
                                reactions.append("rule 6 " + substrate + " " + prod[1:])
                                new_substrates.append(prod[1:])
                    if re.findall('Ab', substrate) == [] and \
                                re.findall('(Ma6Ma6)Mb4', substrate) == []:
                        if re.findall(self.regex_rules[25], substrate) != []:
                            ends = [m.end() for m in re.finditer(self.regex_rules[25], substrate)]
                            for end in ends:
                                prod = substrate[:end] + "(GNb4)" + substrate[end:]
                                reactions.append("rule 10 " + substrate + " " + prod[1:])
                                new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[13], substrate) != []:
                    ends = [m.end() for m in re.finditer(self.regex_rules[13], substrate)]
                    for end in ends:
                        end_split = end - 5
                        prod = substrate[:end_split] + "GNb3" + substrate[end_split:]
                        reactions.append("rule 13 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[4], substrate) != []:
                    inds = [(m.start(),m.end()) for m in re.finditer(self.regex_rules[13], substrate)]
                    for ind in inds:
                        ind_split = ind[0] + 1
                        print substrate[:ind_split], substrate[ind[1]:]
                        prod = substrate[:ind[0]+1] + "NNa3" + substrate[ind[1]:]
                        reactions.append("rule 15 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                #need constraint
                if re.findall(self.regex_rules[16], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[16], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = substrate[:start_split] + "Ab4" + substrate[start_split:]
                        reactions.append("rule 14 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                #same constraint as above
                if re.findall(self.regex_rules[16], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[16], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = substrate[:start_split] + "Ab3" + substrate[start_split:]
                        reactions.append("rule 18 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall('GNb2[\(.*\)|\)]*Ma3', substrate) != [] and \
                                                re.findall(self.regex_rules[2], substrate) == []:
                    if re.findall('[0-9]A[ab]', substrate) == []:
                        starts = [m.start() for m in re.finditer(self.regex_rules[17], substrate)]
                        for start in starts:
                            start_split = start + 4
                            prod = substrate[:start_split] + "(Fa6)" + substrate[start_split:]
                            reactions.append("rule 7 " + substrate + " " + prod[1:])
                            new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[8], substrate) != []:
                    ends = [m.end() for m in re.finditer(self.regex_rules[8], substrate)]
                    for end in ends:
                        end_split = end - 2
                        prod = substrate[:end_split] + "(GNb6)" + substrate[end_split:]
                        reactions.append("rule 16 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[18], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[18], substrate)]
                    for start in starts:
                        start_split = start + 4
                        prod = substrate[:start_split] + "(ANa3)" + substrate[start_split:]
                        prod_other = substrate[:start_split] + "(Aa3)" + substrate[start_split:]
                        reactions.append("rule 25 " + substrate + " " + prod[1:])
                        reactions.append("rule 26 " + substrate + " " + prod_other[1:])
                        new_substrates.append(prod[1:])
                        new_substrates.append(prod_other[1:])
                if re.findall(self.regex_rules[19], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[19], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = substrate[:start_split] + "(Aa3)" + substrate[start_split:]
                        reactions.append("rule 19 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[20], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[20], substrate)]
                    for start in starts:
                        start_split = start + 3
                        prod = substrate[:start_split] + "(Fa4)" + substrate[start_split:]
                        reactions.append("rule 20 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[21], substrate) != []:
                    inds = [(m.start(),m.end()) for m in re.finditer(self.regex_rules[21], substrate)]
                    for ind in inds:
                        start_split = ind[0] + 1
                        prod_before = substrate[:start_split] + "Fa3(" + substrate[start_split:]
                        end_split = ind[1] + 4 - 3
                        prod = prod_before[:end_split] + ")" + prod_before[end_split:]
                        reactions.append("rule 21 " + substrate + " " + prod[1:])
                        reactions.append("rule 24 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[22], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[22], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = substrate[:start_split] + "Fa2" + substrate[start_split:]
                        reactions.append("rule 22 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[23], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[23], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = substrate[:start_split] + "Fa2" + substrate[start_split:]
                        reactions.append("rule 23 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[7], substrate) != []:
                    starts = [m.start() for m in re.finditer(self.regex_rules[7], substrate)]
                    for start in starts:
                        start_split = start + 1
                        prod = substrate[:start_split] + "NNa6" + substrate[start_split:]
                        reactions.append("rule 17 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])
                if re.findall(self.regex_rules[24], substrate) != []:
                    ends = [m.end() for m in re.finditer(self.regex_rules[24], substrate)]
                    for end in ends:
                        end_split = end - 2
                        prod = substrate[:end_split] + "(NNa6)" + substrate[end_split:]
                        reactions.append("rule 27 " + substrate + " " + prod[1:])
                        new_substrates.append(prod[1:])

                count_reactions = 0
                for sub in new_substrates:
                    sub_mass = 0
                    is_capital = False
                    for char in sub:
                        if char.isupper() and not is_capital:
                            sub_mass += 1
                            is_capital = True
                        else:
                            is_capital = False
                    if sub_mass <= self.max_mass:
                        f.write(reactions[count_reactions] + "\n")
                        substrates.append(sub)
                    count_reactions += 1
                count += 1
        f.close()
        return substrates, len(seen_substrates)

mod = glymmer_module([0.5, 0.6], 30)
print "hi"
print mod.generate_network(["Ma2Ma2Ma3(Ma2Ma3(Ma2Ma6)Ma6)Mb4GNb4GN;Asn", "Ma2Ma2Ma3(Ma3(Ma2Ma6)Ma6)Mb4GNb4GN;Asn"])
