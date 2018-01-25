import re

class glymmer_module:

    def __init__(self, parameters, max_mass):
        self.parameters = parameters
        self.max_mass = max_mass
        self.regex_rules = []
        self.regex_rules.append('Ma3\(.*Ma2.*Ma6\)Ma6') #rule 2
        self.regex_rules.append('\(Ma2Ma') #rule 1
        self.regex_rules.append('[ab][0-9]\(GNb[0-9]\)\(.*\)') #rule 5, bisecting glcnac
        self.regex_rules.append('\(Ma3\(Ma3\(Ma6\)Ma6\)Mb4') #rule 8
        self.regex_rules.append('\(Ab4GN') #rule 15
        self.regex_rules.append('\(Ab4GN') #rule 24
        self.regex_rules.append('b4GN') #rule 19
        self.regex_rules.append('\(Ab4GN') #rule 17
        self.regex_rules.append('\(Ab4GNb3A') #rule 16
        self.regex_rules.append('[ab][0-9]\(GNb[0-9]\)\(.*\)') #rule 18

    def generate_network(self, substrates):
        for substrate in substrates:
            substrate = "(" + substrate
            if "Ga3" not in substrate:
                search_string = ""
                #Man1 - specific structure
                if re.findall(self.regex_rules[0], substrate) != []:
                    ind = re.search(self.regex_rules[0], substrate).span()
                    search_string = substrate[ind[0]:ind[1]]
                    print "rule 2"
                #Man1 - general removal of Ma2
                else:
                    if re.findall(self.regex_rules[1], substrate) != []:
                        ind = re.search(self.regex_rules[1], substrate).span()
                        search_string = substrate[ind[0]:ind[1]]
                        print "rule 1"
                # actually removes Ma2
                if search_string != "":
                    str_removal_ind = re.search('Ma2', search_string).span()
                    search_after_removal = search_string[:str_removal_ind[0]] + \
                                            search_string[str_removal_ind[1]:]
                    prod = substrate[:ind[0]] + search_after_removal + substrate[ind[1]:]
                    print "rule 2", prod[1:]
                    substrates.append(prod[1:])
            #GnT1 action
            if re.findall(self.regex_rules[3], substrate) != []:
                ind = re.search(self.regex_rules[3], substrate).span()
                prod = substrate[ind[0]] + "GNb2" + substrate[1:]
                print "rule 8", prod[1:]
                substrates.append(prod[1:])
                # GnT11 action
                # check for bisecting glcnac
                if re.findall(self.regex_rules[2], substrate) == []:
                    ###
            if re.findall(self.regex_rules[4], substrate) != [] or \
                                            substrate.startswith("Ab4GN"):
                #print "rule 15"
                if substrate.startswith("Ab4GN"):
                    prod = "NNa3" + substrate
                    substrates.append(prod)
                else:
                    ind = re.search(self.regex_rules[4], substrate).span()
                    prod = substrate[0] + "NNa3" + substrate[1:]
                    substrates.append(prod)
            if re.findall(self.regex_rules[5], substrate) != [] or \
                                            substrate.startswith("Ab4GN"):
                #print "rule 24" - discrepancy between doc and spreadsheet iirc
            if re.findall(self.regex_rules[6], substrate) != []:
                #print "rule 19" - discrepancy between doc and spreadsheet iirc
            if re.findall(self.regex_rules[7], substrate) != [] or \
                                            substrate.startswith("Ab4GN"):
                #ask about free/outer vs any galactose
                if substrate.startswith("Ab4GN"):
                    prod = "NNa6" + substrate
                    substrates.append(prod)
                else:
                    ind = re.search(self.regex_rules[7], substrate).span()
                    prod = substrate[0] + "NNa6" + substrate[1:]
                    substrates.append(prod)
            if re.findall(self.regex_rules[8], substrate) != []:
                #print "rule 16"
            if re.findall(self.regex_rules[9], substrate) != []:
                #print "rule 18"
                ind = re.search(self.regex_rules[9], substrate).span()
                sub = substrate[:ind[0]] + substrate[ind[1]:]
                if sub.startswith("GN"):
                    substrates.append("Ab3" + substrate)
                if re.findall("\(GN", sub) != []:
                    ind_gn = re.search("\(GN", sub).span()
                    sub = sub[:ind_gn[0]+2] + "Ab3" + sub[ind_gn[0]+2:]
                    #substrates.append(s)
            else:
                #print "rule 18"

        return substrates
