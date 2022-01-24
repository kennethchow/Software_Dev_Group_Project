from website import create_app

# ======= App Creation (from __init__.py ======= #
app = create_app()

# ======= Run File ======= #
if __name__ == '__main__':
    app.run(debug=True)

