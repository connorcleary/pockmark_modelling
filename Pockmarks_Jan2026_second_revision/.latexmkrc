# --- Engine & Output Format Configurations ---
# Force latexmk to output a PDF using pdflatex.
$pdf_mode = 1;

# Specify the base pdflatex compilation command with default flags.
$pdflatex = 'pdflatex -interaction=nonstopmode %O %S';


# --- Directory Routing Configurations ---
# Redirect all auxiliary files (.aux, .log, .toc, etc.) to the 'build' folder.
$aux_dir = 'out-tex';

# Redirect the final output files (.pdf) to the 'build' folder.
$out_dir = 'out-tex';


# --- Bibliography Configurations ---
# Explicitly use BibTeX for bibliography generation.
$bibtex = 'bibtex %O %S';


# --- Glossary & Nomenclature Configurations ---
# Define a custom rule to generate glossaries via the makeglossaries program.
# This tells latexmk to check the input (.glo) and rebuild the output (.gls).
add_cus_dep('glo', 'gls', 0, 'run_makeglossaries');

# Specify the exact command string that executes the makeglossaries rule.
sub run_makeglossaries {
    # Run makeglossaries using the auxiliary directory path so it finds the .glo file.
    return system("makeglossaries -d \"$aux_dir\" \"$_[0]\"");
}


# --- Cleanup Configurations ---
# Ensure custom glossary extension files are deleted during a 'latexmk -c' cleanup.
push @generated_exts, 'glo', 'gls', 'glg', 'ist';