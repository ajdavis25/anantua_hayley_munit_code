#!/usr/bin/env python3
import sys, nltk
from nltk.corpus import wordnet as wn
from wordfreq import zipf_frequency

# make sure NLTK can find downloaded data
nltk.data.path.append('/work/vmo703/nltk_data')

# check command-line argument
if len(sys.argv) < 2:
    print("usage: python lexicon.py <word>")
    sys.exit(1)

word = sys.argv[1].lower()
print(f"\nsearching synonyms for '{word}'...\n")

# collect synonyms with their definitions
synonyms = {}
for syn in wn.synsets(word):
    pos = syn.pos()
    definition = syn.definition()
    for lemma in syn.lemmas():
        name = lemma.name().replace("_", " ")
        if name not in synonyms:
            synonyms[name] = definition

if not synonyms:
    print(f"no WordNet entries found for '{word}'.")
    sys.exit(0)

# rank by rarity (lower frequency = rarer)
ranked = sorted(synonyms.items(), key=lambda kv: zipf_frequency(kv[0], "en"))

print(f"{'word':20} {'ZipfFreq':>8}  definition")
print("-" * 80)
for name, definition in ranked:
    freq = zipf_frequency(name, "en")
    print(f"{name:20} {freq:8.2f}  {definition}")
