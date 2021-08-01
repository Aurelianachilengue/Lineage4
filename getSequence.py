#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse


""" This script extract fasta sequences using IDs in two separated files and save in another file.

#Note: Create more sophisticated logic for matching IDs.""" 


def get_keys(args):
    """This function will save key file into a list."""

    with open(args.keyfile, "r") as kfh:
        keys = [line.rstrip("\n").lstrip(">") for line in kfh]
    return keys


def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Retrieve one or more fastas from a given fasta file."
        )
        parser.add_argument(
            "-f",
            "--fasta",
            action="store",
            required=True,
            help="The multifasta to search.",
        )
        parser.add_argument(
            "-k",
            "--keyfile",
            action="store",
            help="A file of header strings to search the fasta for.",
        )
        parser.add_argument(
            "-s",
            "--string",
            action="store",
            help="Provide a string to look for directly, instead of a file.",
        )
        parser.add_argument(
            "-o",
            "--outfile",
            action="store",
            default=None,
            help="Output file to store the new fasta sequences in. Just prints to screen by default.",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            help="Set whether to print out the key list before the fasta sequences.",
        )
        parser.add_argument(
            "-i",
            "--invert",
            action="store_true",
            help="Invert the search, and retrieve all sequences NOT specified in the keyfile.",
        )
        parser.add_argument(
           "-m",
           "--method",
           action="store",
           choices=["exact", "partial"],
           default="exact",
           help="Search the headers as exact matches, or as partial substring matches. "
       )

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()


def main():
    """Takes a string or list of strings in a text file (one per line) and retreives them and their sequences from a provided fasta file."""
    args = get_args()
    # Call getKeys() to create the list of keys from the provided file:
    if not (args.keyfile or args.string):
        sys.stderr.write("No key source provided. Exiting.")
        sys.exit(1)
    if args.keyfile:
        keys = get_keys(args)
    else:
        keys = args.string.split(",")

    if args.verbose:
        if args.invert is False:
            sys.stderr.write("Fetching the following keys:\n")
            for key in keys:
                sys.stderr.write(key + "\n")
        elif args.invert is True:
            sys.stderr.write(
                "Ignoring the following keys, and retreiving everything else from: {}\n".format(
                    args.fasta
                )
            )
            for key in keys:
                sys.stderr.write(key + "\n")
        sys.stderr.write(
            "-" * 80 + "\n"
        )

    # Parse in the multifasta and assign an iterable variable:
    to_write = []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        if args.invert is False:
            if args.method == "exact":
                if rec.id in keys:
                    print(rec.format("fasta"))
                    to_write.append(rec)

            elif args.method == "partial":
                if any(key in rec.description for key in keys):
                    print(rec.format("fasta"))
                    to_write.append(rec)

        elif args.invert is True:
            if args.method == "exact":
                if rec.id not in keys:
                    print(rec.format("fasta"))
                    to_write.append(rec)
            elif args.method == "partial":
                if all(key not in rec.description for key in keys):
                    print(rec.format("fasta"))
                    to_write.append(rec)

    if args.outfile:
        SeqIO.write(to_write, args.outfile, "fasta")

if __name__ == "__main__":
    main()
