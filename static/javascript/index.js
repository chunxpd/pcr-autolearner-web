var min_seq_num = 1
function add_file(max_file_num) {
    for (let iii = min_seq_num; iii <= max_file_num; iii++) {
        if ($("#in_primer_sequence" + iii).val()) {
            $("#in_primer_sequence" + iii).show();
            $("#in_primer_sequence_label" + iii).show();
            $("#in_primer_sequence_chbox" + iii).show();
        }
        if ($("#in_primer_sequence" + iii).is(":hidden")) {
            $("#in_primer_sequence" + iii).show();
            $("#in_primer_sequence_label" + iii).show();
            $("#in_primer_sequence_chbox" + iii).show();
            break;
        }
    }
};

function del_file(max_file_num) {
    for (let iii = max_file_num; iii >= min_seq_num; iii--) {
        if ($("#in_primer_sequence" + iii).val()) {
            break;
        }
        if ($("#in_primer_sequence" + iii).is(":visible")) {
            $("#in_primer_sequence" + iii).hide();
            $("#in_primer_sequence_label" + iii).hide();
            $("#in_primer_sequence_chbox" + iii).hide();
            break;
        }
    }
};


function rm_whitespace(obj_name) {
    var a = $("#" + obj_name).val().replace(/ /gi, '');
    $("#" + obj_name).val(a);
};

function populate_example() {
    $("#in_primer_sequence1").val("CAGGGCATCGCGCACGCAT")
    $("#in_primer_sequence2").val("TCCACGGCCTGGGCTACCA")
    $("#in_primer_sequence3").val("CACTCCGCCGATCGTTTTTCC")
    $("input[type=checkbox]").prop("checked", false)
    $("#in_primer_sequence2_reverseq").prop("checked", true)
    $("#in_exclude_txids").val("Mycobacterium (taxid: 1763)")
    add_file(100);

};

function populate_example2() {
    $("#in_primer_sequence1").val("CGGCTTCTCCGAGATCATGATG")
    $("#in_primer_sequence2").val("GAAGAACAGGTTCTCCAGCAGG")
    $("#in_primer_sequence3").val("TGCTSGACATCTACCGVAARC")
    $("input[type=checkbox]").prop("checked", false)
    $("#in_primer_sequence1_reverseq").prop("checked", true)
    $("#in_primer_sequence3_reverseq").prop("checked", true)
    $("#in_exclude_txids").val("Mycobacterium (taxid: 1763)")
    add_file(100);
};
