## Getting Started
Clone the repository as follows:
```bash
git clone https://github.com/gioaca00/roche
cd roche
```
Create a conda environment and install the required packages:
```bash
conda env create -n roche-env python=3.11 -y
conda activate roche-env
pip install poetry
poetry install
```

To start using the chatbot, just run the following command:
```bash
python main.py
```

The query is going to be created automatically by the chatbot according to PubMed query construction, following the instructions given by the user. The chatbot will ask to the user to enter a query, followed by the filters that the user wants to apply. The chatbot will then return the 5 articles that best match the query, writing the title, the abstract summarized by an LLM and the first author of the article. \\
The user can also ask for the full abstract of more details about the article.