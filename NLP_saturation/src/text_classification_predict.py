#Run a text classifier. This script is adapted from https://huggingface.co/docs/transformers/main/en/tasks/sequence_classification

import sys
import os
import torch
from datasets import load_dataset
from transformers import pipeline
from transformers import AutoTokenizer
from transformers import AutoModelForSequenceClassification

torch.use_deterministic_algorithms(True)
torch.manual_seed(2025)
model_dir = sys.argv[1]
text_file = sys.argv[2]

def chew_tokens(text, tokenizer, cap = 500):
    segments = [text]
    returnSegments = []
    while segments:
        seg = segments[0]
        segments = segments[1:]
        words = seg.split(' ')
        tokenized_input = tokenizer(words, is_split_into_words=True)
        if len(tokenized_input["input_ids"]) <= cap:
            returnSegments.append(seg)
            continue
        sentences = seg.split(". ")
        if len(sentences) > 1 and not(len(sentences) == 2 and sentences[-1] == ""):
            halfway = int(len(sentences) / 2)
            half1 = ". ".join(sentences[:halfway]) + '. '
            half2 = ". ".join(sentences[halfway:])
        elif len(words) < 2 or (len(words) == 2 and words[-1] == ""):
            halfway = int(len(seg) / 2)
            half1 = seg[:halfway]
            half2 = seg[halfway:]
        else:
            halfway = int(len(words) / 2)
            half1 = ' '.join(words[:halfway]) + ' '
            half2 = ' '.join(words[halfway:])
        segments = [half1, half2] + segments
    return returnSegments


classifier = pipeline("sentiment-analysis", model=model_dir)
for line in open(text_file):
    line = line.strip()
    line = line.split("\t")
    if len(line) < 2: continue
    if not line[0] or not line[1]: continue
    tokenizer = AutoTokenizer.from_pretrained(model_dir)
    model = AutoModelForSequenceClassification.from_pretrained(model_dir)
    model.eval()
    for segment in chew_tokens(line[1], tokenizer):
        inputs = tokenizer(segment, return_tensors="pt")
        with torch.no_grad(): logits = model(**inputs).logits
        logits = list(logits[0])
        for i in range(len(logits)): logits[i] = float(logits[i])
        predicted_class_id = 0
        champval = logits[0]
        for i in range(1, len(logits)):
            if logits[i] > champval:
                predicted_class_id = i
                champval = logits[i]
        for i in range(len(logits)): logits[i] = str(logits[i])
        logits = "\t".join(logits)
        print(line[0] + "\t" + model.config.id2label[predicted_class_id] + "\t" + logits + "\t" + line[1], flush=True)

