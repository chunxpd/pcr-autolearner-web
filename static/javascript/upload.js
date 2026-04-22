// for uploading fasta files
var uploaded_files;
$(document).ready(function () {
    $('#uploadFastas').submit(function (event) {
        $("#loader-icon").show();
        if ($('#uploadFile1').val()) {
            event.preventDefault();
            $(this).ajaxSubmit({
                target: '#uploadResults',
                async: true,
                beforeSubmit: function () {
                    $('.progress-bar').width('0%');
                },
                uploadProgress: function (event, position, total, percentageComplete) {
                    $('.progress-bar').animate({
                        width: percentageComplete + '%'
                    }, {
                        duration: 0
                    });
                },
                success: function (data) {
                    $('.progress-bar').width('0%');
                    $('#uploadResults').show();
                    var res = "<h4>Uploaded files:</h4><ul>"
                    for (let iii = 0; iii < data["size"]; iii++) {
                        res += "<li>" + data["data"][iii] + "</li>"
                    }
                    res += "</ul>"
                    $('#uploadResults').html(res);
                    $('#uploadResults').append("jobid: " + data["jobid"]);
                    uploaded_files = data;
                },
                resetForm: true
            });
        }
        return false;
    });
});

function add_file(max_file_num) {
    for (let iii = 1; iii <= max_file_num; iii++) {
        if ($("#uploadFile" + iii).is(":hidden")) {
            $("#uploadFile" + iii).show();
            break;
        }
    }
};

function del_file(max_file_num) {
    for (let iii = max_file_num; iii >= 1; iii--) {
        if ($("#uploadFile" + iii).is(":visible")) {
            $("#uploadFile" + iii).hide();
            break;
        }
    }
};
