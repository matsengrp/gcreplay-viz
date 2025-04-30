import os, sys
import subprocess
import argparse
import re
from pathlib import Path
from collections import defaultdict
from pprint import pp


### arg parser ###


class Parser:
    def __init__(self, arg_parser, arg_help={}, arg_default={}, arg_options={}):
        self.build(arg_parser, arg_help, arg_default, arg_options)

    def build(self, arg_parser, arg_help={}, arg_default={}, arg_options={}):
        self.arg_parser = arg_parser
        self.arg_help = arg_help
        self.arg_default = arg_default
        self.arg_options = arg_options

        self.arg_parser.set_defaults(**self.arg_default)
        for action in self.arg_parser._actions:
            if action.help is None:
                action.help = ""
            if action.dest in self.arg_help.keys():
                action.help += f"{self.arg_help[action.dest]} "
            elif action.dest in self.arg_default.keys():
                action.help += f"(default: '{self.arg_default[action.dest]}') "
        return self.arg_parser

    def parse_args(self, commandline_args=sys.argv[1:]):
        if isinstance(commandline_args, str):
            self.commandline_args = commandline_args.split()
        assert isinstance(commandline_args, list)
        self.commandline_args = commandline_args

        args = self.arg_parser.parse_args(self.commandline_args)
        args = defaultdict(lambda: None, vars(args))
        return args

    @staticmethod
    def parse_list(parse_fn=str):
        def parser(args):
            output = [parse_fn(arg) for arg in args.split(",")]
            return output
        return parser

    @staticmethod
    def parse_option(options, parse_fn=str):
        def parser(arg):
            output = parse_fn(arg)
            if output not in options:
                raise argparse.ArgumentTypeError(f"Invalid option: '{arg}' (expected one of {options})")
            return output
        return parser

    @staticmethod
    def parse_option_list(options, parse_fn=str):
        return Parser.parse_list(Parser.parse_option(options, parse_fn))

    @staticmethod
    def parse_input_file():
        def parser(arg):
            path = Path(arg)
            if not path.is_file():
                raise argparse.ArgumentTypeError(f"Input file does not exist (or is not a file): {path}")
            return str(path)
        return parser

    @staticmethod
    def parse_output_file():
        def parser(arg):
            path = Path(arg)
            if not path.parent.exists():
                raise argparse.ArgumentTypeError(f"Output file's parent directory does not exist: {path.parent}")
            return str(path)
        return parser

    @staticmethod
    def parse_input_dir():
        def parser(arg):
            path = Path(arg)
            if not path.is_dir():
                raise argparse.ArgumentTypeError(f"Input directory does not exist (or is not a directory): {path}")
            return str(path)
        return parser

    @staticmethod
    def parse_output_dir():
        def parser(arg):
            path = Path(arg)
            if not path.parent.exists():
                raise argparse.ArgumentTypeError(f"Output directory's parent does not exist: {path.parent}")
            return str(path)
        return parser

    @staticmethod
    def parse_flag(default=True):
        def parser(arg=None):
            if isinstance(arg, bool):
                return arg
            arg = arg.lower()
            if arg in {"t", "true", "y", "yes", "1"}:
                return True
            elif arg in {"f", "false", "n", "no", "0"}:
                return False
            raise argparse.ArgumentTypeError(f"Invalid boolean arg: '{arg}'. Expected format: true/false")
        return parser

    @staticmethod
    def parse_colorhex():
        def parser(arg):
            arg = arg.upper()
            match = re.fullmatch(r'#?([0-9A-Fa-f]{6})', arg)
            if not match:
                raise argparse.ArgumentTypeError(f"Invalid color format: '{arg}'. Expected format: #RRGGBB or RRGGBB.")
            return f"#{match.group(1).upper()}"
        return parser

    @staticmethod
    def parse_range():
        def parser(args):
            output = []
            args = args.split(",")
            for arg in args:
                arg = arg.split("-")
                print(f"arg: len:{len(arg)} {arg} ")
                if len(arg) == 1:
                    output.append([arg])
                elif len(arg) == 2:
                    output.append([arg[0], arg[1]])
                else:
                    raise argparse.ArgumentTypeError(f"Invalid range format: '{arg}'. Expected format: 1,4-5,6,10")
            return output
        return parser

    @staticmethod
    def flatten_range(int_ranges, inclusive=True):
        output = []
        for int_pair in int_ranges:
            if len(int_pair) == 1:
                output.append(int_pair)
            elif len(int_pair) == 2:
                output += range(int_pair[0], int_pair[1] + inclusive)
        return output


### utilities ###


def run_command(command, do_print=True):
    if do_print:
        print(f"COMMAND: {command}")
    try:
        output = subprocess.run(
            command, shell=True, check=True, capture_output=True, text=True
        )
    except Exception as e:
        print(f"COMMAND failed with exception: {e}")
        return None
    else:
        print(f"COMMAND successful!")
        return output


class ColorPrinter:
    class colors:
        BLACK = "\033[30m"
        RED = "\033[91m"
        GREEN = "\033[92m"
        YELLOW = "\033[93m"
        BLUE = "\033[94m"
        MAGENTA = "\033[95m"
        CYAN = "\033[96m"
        WHITE = "\033[97m"

        BG_BLACK = "\033[40m"
        BG_RED = "\033[41m"
        BG_GREEN = "\033[42m"
        BG_YELLOW = "\033[43m"
        BG_BLUE = "\033[44m"
        BG_MAGENTA = "\033[45m"
        BG_CYAN = "\033[46m"
        BG_WHITE = "\033[47m"

        BOLD = "\033[1m"
        UNDERLINE = "\033[4m"
        RESET = "\033[0m"

    @staticmethod
    def bash_color_code(R=255, G=255, B=255):
        code = f"\033[38;2;{R};{G};{B};0m"
        return code

    @staticmethod
    def print(*args, color=None, bg_color=None, style=None, sep=" ", end="\n"):
        color_code = ""
        if color:
            color_code += color
        if bg_color:
            color_code += bg_color
        if style:
            color_code += style
        color_reset = ""
        if color_code != "":
            color_reset = ColorPrinter.colors.RESET
        message = sep.join(map(str, args))
        print(f"{color_code}{message}{color_reset}", end=end)

    @staticmethod
    def set_color(color):
        print(f"{color}")

    @staticmethod
    def unset_color():
        print(f"{ColorPrinter.colors.RESET}")


colors = ColorPrinter.colors
cprint = ColorPrinter.print
cprint_set_color = ColorPrinter.set_color
cprint_unset_color = ColorPrinter.unset_color


class Encoder:
    long2short_dict = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }

    short2long_dict = {v: k for k, v in long2short_dict.items()}

    @staticmethod
    def long2short(long_code):
        if long_code in Encoder.long2short_dict.keys():
            return Encoder.long2short_dict[long_code]
        return "X"

    @staticmethod
    def short2long(short_code):
        return Encoder.short2long_dict[short_code]


### unit tests ###


def build_test_parser():
    arg_parser = argparse.ArgumentParser("test parser")
    arg_parser.add_argument("--run", type=Parser.parse_option_list(['build', 'join']))
    arg_parser.add_argument("--input-dir", type=Parser.parse_input_dir())
    arg_parser.add_argument("--output-dir", type=Parser.parse_output_dir())
    arg_parser.add_argument("--input-file", type=Parser.parse_input_file())
    arg_parser.add_argument("--output-file", type=Parser.parse_output_file())
    arg_parser.add_argument("--names", type=Parser.parse_list(parse_fn=str))
    arg_parser.add_argument("--numbers", type=Parser.parse_list(parse_fn=int))
    arg_parser.add_argument("--color", type=Parser.parse_colorhex())
    arg_parser.add_argument("--colors", type=Parser.parse_list(Parser.parse_colorhex()))
    arg_parser.add_argument("--range", type=Parser.parse_range())
    arg_parser.add_argument("--flag", type=Parser.parse_flag(), const=True, default=False)

    parser = Parser(arg_parser)
    return parser


def main():
    parser = build_test_parser()
    args = parser.parse()

    print(f'args:')
    pp(dict(args))


if __name__ == "__main__":
    main()
