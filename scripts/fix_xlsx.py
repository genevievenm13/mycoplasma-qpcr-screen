"""Repair Bio-Rad CFX .xlsx exports (backslash paths, lowercase content-types) so pandas can read them."""
import sys, zipfile, pathlib
src, dst = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
dst.mkdir(parents=True, exist_ok=True)
for f in sorted(src.glob("*.xlsx")):
    with zipfile.ZipFile(f) as zin, zipfile.ZipFile(dst / f.name, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            name = info.filename.replace("\\", "/")
            if name.lower() == "[content_types].xml":
                name = "[Content_Types].xml"
            if name.lower() == "xl/sharedstrings.xml":
                name = "xl/sharedStrings.xml"
            zout.writestr(name, zin.read(info.filename))
    print("fixed", f.name)
