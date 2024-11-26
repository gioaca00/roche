# Load model directly
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from torch import cuda


def load_model():
    tokenizer = AutoTokenizer.from_pretrained("Falconsai/medical_summarization")
    model = AutoModelForSeq2SeqLM.from_pretrained("Falconsai/medical_summarization")
    return tokenizer, model


def summarize_abstract(text, min_length=32, max_length=256):
    """Summarize a given text using a pre-trained model."""
    device = "cuda" if cuda.is_available() else "cpu"
    tokenizer, model = load_model()
    model = model.to(device=device)
    inputs = tokenizer(text, return_tensors="pt", max_length=1024, truncation=True).to(
        model.device
    )
    outputs = model.generate(
        **inputs, max_length=max_length, min_length=min_length, num_return_sequences=1
    )
    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    return generated_text
