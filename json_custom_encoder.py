"""=================================================================================================
Format for json output
Modified from https://stackoverflow.com/questions/16264515/json-dumps-custom-formatting

Michael Gribskov     21 February 2024
================================================================================================="""
import json

class CompactJSONEncoder(json.JSONEncoder):
    """A JSON Encoder that puts small lists on single lines."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.indentation_level = 0

    def encode(self, o):
        """Encode JSON object *o* with respect to single line lists."""

        if isinstance(o, (list, tuple)):
            if self._is_single_line_list(o):
                return "[" + ", ".join(json.dumps(el) for el in o) + "]"
            else:
                self.indentation_level += 1
                output = [self.indent_str + self.encode(el) for el in o]
                self.indentation_level -= 1
                return "[\n" + ",\n".join(output) + "\n" + self.indent_str + "]"

        elif isinstance(o, dict):
            self.indentation_level += 1
            output = [self.indent_str + f"{json.dumps(k)}: {self.encode(v)}" for k, v in o.items()]
            self.indentation_level -= 1
            return "{\n" + ",\n".join(output) + "\n" + self.indent_str + "}"

        else:
            return json.dumps(o)

    def _is_single_line_list(self, o):
        """-----------------------------------------------------------------------------------------
        This controls what should be a single line

        :param o: list      list to format
        :return:
        -----------------------------------------------------------------------------------------"""
        n_items = 5
        line_len = 120
        if isinstance(o, (list, tuple)):
            return not any(isinstance(el, (list, tuple, dict)) for el in o) \
                and len(o) <= n_items \
                and len(str(o)) - 2 <= line_len

    @property
    def indent_str(self) -> str:
        return " " * self.indentation_level * self.indent

    def iterencode(self, o, **kwargs):
        """Required to also work with `json.dump`."""
        return self.encode(o)
# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    exit(0)
