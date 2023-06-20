import re
import argparse

"""
Remove all instances of Amyloid-Beta strcutures from the MSA 
step of AlphaFold so we can oberserve the results in isolation
"""

parser = argparse.ArgumentParser(
    description="Removes all instances of protiens that match a regex string from a HMMER database query"
)
parser.add_argument(
    "-i",
    "--input",
    nargs=1,
    help="Input file (typically a .sto file such as pdb_hits.sto)",
    required=True,
)
parser.add_argument("-o", "--output", nargs=1, help="Output file", required=True)
parser.add_argument(
    "--regex",
    nargs=1,
    default=[r"a.*beta|amyloid"],
    help="Optional regex for parsing (defaults to removing AB results)",
)
args = parser.parse_args()


class HMMERReader:
    def __init__(self, fname: str):
        self.datafile = [l.strip() for l in open(fname, "r").readlines()]
        self._contents = None
        self.init_size = (len(self.datafile) - 5) // 3
        self.pr = None

        self.header = []
        self.body = []
        self.tail = []

        pointer = 0
        # breaks between header, body, and tail are seperated by blank likes
        while self.datafile[pointer]:
            self.header.append(self.datafile[pointer])
            pointer += 1
        
        pointer += 1
        while self.datafile[pointer]:
            self.body.append(self.datafile[pointer])
            pointer += 1
        
        pointer += 1
        self.tail = self.datafile[pointer:len(self.datafile)]

        self.pruned_body = []
        self.pruned_tail = []

    # Remove the database hits that match given regex
    def prune(self, prune_regex: str):
        prune_regex = prune_regex.lower()
        self.pr = prune_regex
        prots_to_remove = set()

        for line in self.body:
            data = line.split(" ")
            prot_code = data[1]

            # determine proteins whose title match our query and need to be removed
            if re.findall(prune_regex, line.lower()):
                prots_to_remove.add(prot_code)
            # if there is no match, write to file
            else:
                self.pruned_body.append(line)
        
        for line in self.tail:
            # match to GR
            prot_code = None

            if line.startswith("#=GR"):
                data = line.split(" ")
                prot_code = data[1]

            # match to raw title
            elif re.findall(r".*/.*", line.lower()):
                data = line.split(" ")
                prot_code = data[0]

            # if prot code is not one that is matched in our query, add it to our pruned file
            if not prot_code or prot_code not in prots_to_remove:
                self.pruned_tail.append(line)

    # write result to file
    def output(self, fname: str):
        # concat file contents
        lines = self.header + [""] + self.pruned_body + [""] + self.pruned_tail
        print(f"Pruned {len(self.body) - len(self.pruned_body)} Entries with regex r'{self.pr}'")

        with open(fname, "w") as f:
            for line in lines:
                f.write(f"{line}\n")


def main():
    hreader = HMMERReader(args.input[0])
    hreader.prune(args.regex[0])
    hreader.output(args.output[0])


if __name__ == "__main__":
    main()
