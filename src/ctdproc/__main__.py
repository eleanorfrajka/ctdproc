import sys

from .main import CTD_process


def main():

    if len(sys.argv) != 2:
        print(
            "Usage: python -m ctdproc config.yaml"
        )
        sys.exit(1)

    yaml_file = sys.argv[1]

    CTD_process(yaml_file)


if __name__ == "__main__":
    main()
