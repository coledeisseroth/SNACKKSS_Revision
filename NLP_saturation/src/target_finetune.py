#This code is adapted from the tutorial code at https://huggingface.co/docs/transformers/main/en/tasks/token_classification
import sys
import os
from collections import defaultdict
from datasets import load_dataset
from datasets import Dataset
from datasets import DatasetDict
import torch
import numpy as np

np.random.seed(2025)
torch.manual_seed(2025)
torch.use_deterministic_algorithms(True)
train_json = sys.argv[1]
out_dir = sys.argv[2]
ckpt = sys.argv[3]
steps = int(sys.argv[4])
epochs = int(sys.argv[5])
batchsize = int(sys.argv[6])

label_list = ["O"]
for line in ["GENE", "CHEM"]:
    label_list.append("B-" + line)
    label_list.append("I-" + line)

id2label = {}
label2id = {}
for i in range(len(label_list)):
    id2label[i] = label_list[i]
    label2id[label_list[i]] = i


my_dataset = load_dataset('json', data_files={'train': train_json})

from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained(ckpt, add_prefix_space=True)


example = my_dataset["train"][0]
tokenized_input = tokenizer(example["tokens"], is_split_into_words=True)
tokens = tokenizer.convert_ids_to_tokens(tokenized_input["input_ids"])

def tokenize_and_align_labels(examples):
    tokenized_inputs = tokenizer(examples["tokens"], truncation=True, max_length=500, is_split_into_words=True)
    labels = []
    for i, label in enumerate(examples[f"ner_tags"]):
        word_ids = tokenized_inputs.word_ids(batch_index=i)
        previous_word_idx = None
        label_ids = []
        for word_idx in word_ids:
            if word_idx is None:
                label_ids.append(-100)
            elif word_idx != previous_word_idx:
                label_ids.append(label[word_idx])
            elif word_idx == previous_word_idx:
                if id2label[label[word_idx]][0] == "B": label_ids.append(label[word_idx] + 1)
                else: label_ids.append(label[word_idx])
            else:
                label_ids.append(-100)
            previous_word_idx = word_idx
        labels.append(label_ids)

    tokenized_inputs["labels"] = labels
    return tokenized_inputs

tokenized_dataset = my_dataset.map(tokenize_and_align_labels, batched=True)

from transformers import DataCollatorForTokenClassification

data_collator = DataCollatorForTokenClassification(tokenizer=tokenizer)

import evaluate

seqeval = evaluate.load("seqeval")

labels = [label_list[i] for i in example[f"ner_tags"]]

def compute_metrics(p):
    predictions, labels = p
    predictions = np.argmax(predictions, axis=2)

    true_predictions = [
        [label_list[p] for (p, l) in zip(prediction, label) if l != -100]
        for prediction, label in zip(predictions, labels)
    ]
    true_labels = [
        [label_list[l] for (p, l) in zip(prediction, label) if l != -100]
        for prediction, label in zip(predictions, labels)
    ]

    results = seqeval.compute(predictions=true_predictions, references=true_labels)
    return {
        "precision": results["overall_precision"],
        "recall": results["overall_recall"],
        "f1": results["overall_f1"],
        "accuracy": results["overall_accuracy"],
    }

from transformers import AutoModelForTokenClassification, TrainingArguments, Trainer

model = AutoModelForTokenClassification.from_pretrained(ckpt, num_labels=len(label_list), id2label=id2label, label2id=label2id, ignore_mismatched_sizes=True)
model.eval()

training_args = TrainingArguments(
    output_dir=out_dir,
    learning_rate=2e-5,
    per_device_train_batch_size=batchsize,
    num_train_epochs=epochs,
    weight_decay=0.01,
    evaluation_strategy="no",
    #save_strategy="no",
    save_steps = steps,
    seed=2025,
    use_cpu=True,
    gradient_checkpointing=True,
    #load_best_model_at_end=True,
    push_to_hub=False,
    report_to="none"
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset["train"],
    tokenizer=tokenizer,
    data_collator=data_collator,
    compute_metrics=compute_metrics,
)

trainer.train()
trainer.save_model()

