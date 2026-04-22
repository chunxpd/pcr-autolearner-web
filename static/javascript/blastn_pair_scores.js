var isRun = 0;
var cd_hit_files;
// uploaded_files 의 내용
// {
//  "jobid": jobid (directory name where xml files are stored),
//  "size": # of xml files,
//  "data": list of xml files
// }
$(document).ready(function () {
    // for cd-hit : reduction of redundancy
    $("#run_hsp_pair_scores").click(function () {
        if (isRun > 0) {
            // ajax 실행 중에 중복해서 실행 되지 않게 한다.
            alert("blastn HSPs pair score 계산 실행 중");
            return;
        }
        $.ajax({
            type: 'POST',
            // url: '{{url_for("cal_blastn_xml_pair")}}',
            url: '/cal_blastn_xml_pair',
            data: JSON.stringify({
                "files": uploaded_files["data"],
                "jobid": uploaded_files["jobid"],
                "size": uploaded_files["size"]
            }),
            dataType: 'JSON',
            contentType: "application/json",
            success: function (data) {
                // $("#hsp_pair_scores_run_res").html(data["res"]);
                cd_hit_files = data["res"]
                $("#hsp_pair_scores_run_res").html(
                    "<div class='card card-body'>" +
                    "tsv file: " +
                    "<a href=" +
                    "/local_download_general/__JOBID__/__FILENAME__".replace("__JOBID__", data["jobid"]).replace("__FILENAME__", data["filename"]) +
                    '>' +
                    data['filename'] + "</a>" +
                    "</div>");
                // $("#msa_input").html(cd_hit_files);
            },
            beforeSend: function () {
                isRun = isRun + 1; // ajax 실행 중에 중복해서 실행 되지 않게 한다.
            },
            complete: function () {
                isRun = isRun - 1;
            },
            error: function (request, status, error) {
                alert('run hsp pair score error:[' + status + '][' + error + ']');
            }
        })
    });
});
