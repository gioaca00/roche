import questionary

def multiple_choice(options, question):
    """Prompt the user to choose which additional information to provide."""
    selected_option = questionary.select(
        question,
        choices=options,
    ).ask()
    return selected_option
