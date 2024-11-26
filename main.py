from Bio import Entrez
from roche.utils import multiple_choice
from roche.language_model import summarize_abstract


def search_and_fetch(query, max_results=5):
    """Search PubMed for a query and fetch metadata for the top results."""
    with Entrez.esearch(db="pubmed", term=query, retmax=max_results) as handle:
        results = Entrez.read(handle)
    article_ids = results["IdList"]

    try:
        with Entrez.efetch(
            db="pubmed", id=article_ids, rettype="abstract", retmode="xml"
        ) as handle:
            records = Entrez.read(handle)
    except:
        print("Your search did not return any results. Please try again.")
        return

    for i, article in enumerate(records["PubmedArticle"]):
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        first_author = article["MedlineCitation"]["Article"]["AuthorList"][0].get(
            "LastName", ""
        )
        abstract = (
            article["MedlineCitation"]["Article"]
            .get("Abstract", {})
            .get("AbstractText", ["No abstract"])[0]
        )
        if abstract == "":
            summarized_abstract = "No abstract available."
        else:
            summarized_abstract = summarize_abstract(abstract) + "..."
        print(
            f"{i+1}. Title: {title}\n   Authors: {first_author} et al.\n   Abstract: {summarized_abstract}\n"
        )

    return records


def login_and_create_query():
    """Login to PubMed and create a search query based on user input."""
    if not Entrez.email:
        print("Please provide your email for identification.")
        Entrez.email = input("Email: ")
    query = input("Enter your query (e.g., 'COVID-19 vaccine efficacy'): ")
    interface_options = ["No, continue", "Keywords", "Time Range", "Other"]
    interface_question = "Do you want to provide additional information for the search?"
    while True:
        choice = multiple_choice(interface_options, interface_question)
        if choice == "Keywords":
            keywords = input("Enter additional keywords: (e.g., 'cancer treatment') ")
            query += f' AND  ("{keywords}")'
        elif choice == "Time Range":
            date_range = input("Enter time range (e.g., '2021-2023'): ").split("-")
            if len(date_range) == 1:
                query += ' AND ("' + date_range[0] + '"[Date - Publication])'
            else:
                query += (
                    ' AND ("'
                    + date_range[0]
                    + '"[Date - Publication] : "'
                    + date_range[1]
                    + '"[Date - Publication])'
                )
        elif choice == "Other":
            field = input("Enter field name: (e.g., 'Title') ")
            value = input("Enter field value: (e.g., 'COVID-19') ")
            query += f' AND  ("{value}"[{field}])'
        else:
            break
    return query


def more_details(articles):
    """Provide additional details about a specific article."""
    while True:
        article_number = input(
            "Enter the number of the article you want more details about (e.g., 1), or press enter to stop: "
        )
        if article_number == "":
            break
        try:
            article_number = int(article_number)
        except ValueError:
            print("Invalid article number. Please try again.")
        if 1 <= article_number <= len(articles["PubmedArticle"]):
            article_options = [
                "I'm done with this paper",
                "Full Abstract",
                "Authors List",
                "Article Date",
                "Journal Information",
            ]
            article_question = "What do you want to know?"
            while True:
                choice = multiple_choice(article_options, article_question)
                if choice == "Full Abstract":
                    print(
                        articles["PubmedArticle"][article_number - 1][
                            "MedlineCitation"
                        ]["Article"]
                        .get("Abstract", {})
                        .get("AbstractText", ["No abstract"])[0]
                    )
                elif choice == "Authors List":
                    authors = articles["PubmedArticle"][article_number - 1][
                        "MedlineCitation"
                    ]["Article"]["AuthorList"]
                    for author in authors:
                        print(
                            f"{author.get('ForeName', '')} {author.get('LastName', '')}"
                        )
                elif choice == "Article Date":
                    article_date = articles["PubmedArticle"][article_number - 1][
                        "MedlineCitation"
                    ]["Article"]["ArticleDate"]
                    if article_date:
                        print(
                            f"{article_date[0].get('Day', '')}/{article_date[0].get('Month', '')}/{article_date[0].get('Year', '')}"
                        )
                    else:
                        print("No specific date available.")
                elif choice == "Journal Information":
                    journal = articles["PubmedArticle"][article_number - 1][
                        "MedlineCitation"
                    ]["Article"]["Journal"]
                    print(f"Journal Name: {journal.get('Title', '')}")
                else:
                    break
        else:
            print("Invalid article number. Please try again.")


if __name__ == "__main__":
    print("Welcome to PubMed Chatbot!")
    query = login_and_create_query()

    print("Searching PubMed...")
    articles = search_and_fetch(query, max_results=5)

    if not articles:
        exit()

    print("Do you want more details about a specific article?")
    more_details(articles)
