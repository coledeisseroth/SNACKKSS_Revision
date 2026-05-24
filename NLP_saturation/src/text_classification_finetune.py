#This code is adapted from the tutorial code at https://huggingface.co/docs/transformers/main/en/tasks/sequence_classification
import sys
import os
import torch
from datasets import load_dataset
from transformers import AutoTokenizer
from transformers import DataCollatorWithPadding
import evaluate
import numpy as np
from transformers import AutoModelForSequenceClassification, TrainingArguments, Trainer

np.random.seed(2025)
torch.use_deterministic_algorithms(True)
torch.manual_seed(2025)
train_json = sys.argv[1]
out_dir = sys.argv[2]
ckpt = sys.argv[3]
steps = int(sys.argv[4])

my_dataset = load_dataset('json', data_files={'train':train_json})

tokenizer = AutoTokenizer.from_pretrained(ckpt)

def preprocess_function(examples):
    return tokenizer(examples["text"], truncation=True, max_length=500)

tokenized_dataset = my_dataset.map(preprocess_function, batched=True)

data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

accuracy = evaluate.load("accuracy")

def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    predictions = np.argmax(predictions, axis=1)
    return accuracy.compute(predictions=predictions, references=labels)

id2label = {0: "NEGATIVE", 1: "POSITIVE"}
label2id = {"NEGATIVE": 0, "POSITIVE": 1}

model = AutoModelForSequenceClassification.from_pretrained(ckpt, num_labels=2, id2label=id2label, label2id=label2id, ignore_mismatched_sizes=True)
model.eval()

training_args = TrainingArguments(
    output_dir=out_dir,
    learning_rate=2e-5,
    per_device_train_batch_size=1,
    num_train_epochs=2,
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


