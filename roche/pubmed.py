from Bio import Entrez
import questionary

# Provide your email for identification
Entrez.email = "example@example.com"

def multiple_choice(options, question):
    """Prompt the user to choose which additional information to provide."""
    selected_option = questionary.select(
        question,
        choices=options,
        default=options[0]
    ).ask()
    return selected_option

def search_pubmed(query, max_results=5):
    """Search PubMed and return a list of article IDs."""
    with Entrez.esearch(db="pubmed", term=query, retmax=max_results) as handle:
        results = Entrez.read(handle)
    return results["IdList"], results["Count"]

def fetch_details(article_ids):
    """Fetch details of articles given their IDs."""
    with Entrez.efetch(db="pubmed", id=article_ids, rettype="abstract", retmode="xml") as handle:
        records = Entrez.read(handle)
    return records

def search_and_fetch(query, max_results=5):
    """Search PubMed for a query and fetch metadata for the top results."""
    article_ids, _ = search_pubmed(query, max_results)
    return fetch_details(article_ids)

if __name__ == "__main__":
    # interface
    print("Welcome to PubMed Chatbot!")
    query = input("Enter your query (e.g., 'COVID-19 vaccine efficacy'): ")
    interface_options = ["No, continue", "Keywords", "Time Range", "Other"]
    interface_question = "Do you want to provide additional information for the search?"
    while True:
        choice = multiple_choice(interface_options, interface_question)
        if choice == "Keywords":
            query += " " + input("Enter additional keywords: (e.g., 'cancer treatment') ")
        elif choice == "Time Range":
            query += " AND " + input("Enter time range (e.g., '2021-2023'): ")
        elif choice == "Other":
            query += " " + input("Enter additional information: (e.g., author, name or journal) ")
        else:
            break

    # search for articles
    print("Searching PubMed...")
    articles = search_and_fetch(query, max_results=5)
    abstract_length = 200
    for i, article in enumerate(articles["PubmedArticle"]):
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        first_author = article["MedlineCitation"]["Article"]["AuthorList"][0].get("LastName", "")
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No abstract"])[0][:abstract_length]
        print(f"{i+1}. Title: {title}\n   Authors: {first_author} et al.\n   Abstract: {abstract}...\n")
    # let users ask more details about one of the articles
    while True:
        article_number = input("Enter the number of the article you want more details about (e.g., 1), or press enter to stop: ")
        if article_number == "":
            break
        try:
            article_number = int(article_number)
        except ValueError:
            print("Invalid article number. Please try again.")
        if 1 <= article_number <= len(articles["PubmedArticle"]):
            article_options = ["I'm done with this paper", "Full Abstract", "Authors List", "ArticleDate", "Journal Information", "All available info (not formatted)"]
            article_question = "What do you want to know?"
            while True:
                choice = multiple_choice(article_options, article_question)
                if choice == "Full Abstract":
                    print(articles["PubmedArticle"][article_number - 1]["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No abstract"])[0])
                elif choice == "Authors List":
                    authors = articles["PubmedArticle"][article_number - 1]["MedlineCitation"]["Article"]["AuthorList"]
                    for author in authors:
                        print(f"{author.get('ForeName', '')} {author.get('LastName', '')}")
                elif choice == "ArticleDate":
                    article_date = articles["PubmedArticle"][article_number - 1]["MedlineCitation"]["Article"]["ArticleDate"]
                    print(f"{article_date[0].get('Day', '')}/{article_date[0].get('Month', '')}/{article_date[0].get('Year', '')}")
                elif choice == "Journal Information":
                    journal = articles["PubmedArticle"][article_number - 1]["MedlineCitation"]["Article"]["Journal"]
                    print(f"Journal Name: {journal.get('Title', '')}")
                elif choice == "All available info (not formatted)":
                    print(articles["PubmedArticle"][article_number - 1])
                else:
                    break
        else:
            print("Invalid article number. Please try again.")