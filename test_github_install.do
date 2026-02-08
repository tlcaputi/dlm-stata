/*
    test_github_install.do
    ======================
    Tests installing dlm from the GitHub repo via net install.
    This is what real users will do.
*/

clear all
set more off

// ============================================================================
// Uninstall any existing version
// ============================================================================
capture ado uninstall dlm
display as text "Cleaned any previous install"

// ============================================================================
// Install from GitHub
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "Installing dlm from GitHub..."
display as text "{hline 72}"

net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/") replace

// ============================================================================
// Verify it works
// ============================================================================
display as text ""
display as text "{hline 72}"
display as text "Verifying installation..."
display as text "{hline 72}"

which dlm
which dlm_gen_data

// Generate data and run
dlm_gen_data, n_groups(200) n_times(15) treat_prob(0.3) seed(42) clear
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)

display as text ""
display as result "GitHub install test PASSED"
display as text "Users can install with:"
display as text `"  net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/")"'
