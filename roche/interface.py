from pymed import PubMed

def search_pubmed(keywords, mindate=None, maxdate=None, retmax=5):
    """Search PubMed for articles based on keywords and date range."""
    pubmed = PubMed(tool="PubMedSearcher", email="example@example.com")

    keywords = " ".join(keywords)

    if mindate is not None and maxdate is not None:
        date = f"AND ({mindate}:[Date - Create] {maxdate}:[Date - Create])"
    elif mindate is not None:
        date = f"AND ({mindate}:[Date - Create])"
    elif maxdate is not None:
        date = f"AND ({maxdate}:[Date - Create])"
    else:
        date = ""
    query = f"{keywords} {date}"
    articles = pubmed.query(query, max_results=retmax)
    return articles

def process_input(user_input):
    """Process user input and extract query parameters."""
    keywords = user_input.split()
    print("Extracted Keywords:", keywords)
    return {"keywords": keywords, "mindate": None, "maxdate": None}

def display_results(results, abstract_length=150):
    """Display a list of articles in a user-friendly format."""
    print("\nSearch Results:")
    for i, article in enumerate(results, 1):
        print(f"{i}. Title: {article['title']}\n   Authors: {article['authors']}")
        print(f"   Abstract: {article['abstract'][:abstract_length]}...\n")

def main():
    print("Welcome to PubMed Chatbot!")
    user_input = input("Enter your query (e.g., 'Find articles about \"COVID-19 vaccine efficacy\"'): ")
    user_keywords = input("Optionally add keywords (e.g., 'cancer treatment'): ")
    user_timerange = input("Optionally add timerange (e.g., '2021-2023'): ")
    user_filters = input("Optionally add more filters (e.g., author, name or journal).: ")
    
    user_keywords, user_mindate, user_maxdate = process_input(user_input)

    print("Searching PubMed...")
    articles = search_pubmed(user_keywords, retmax=5)   

    sample_results = []
    for article in articles:
        sample_results.append({
            "title": article.title,
            "authors": f"{article.authors[0]['initials']} {article.authors[0]['lastname']} et al.",
            "abstract": article.abstract,
        })
    
    display_results(sample_results, 150)

if __name__ == "__main__":
    main()
