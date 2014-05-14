{%- extends 'html_basic.tpl' -%}

{% block input_group -%}
<div class="input_hidden">
{{ super() }}
</div>
{% endblock input_group %}

{%- block header -%}
{{ super() }}

<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js"></script>

<style type="text/css">
//div.output_wrapper {
//  margin-top: 0px;
//}
.input_hidden {
  display: none;
//  margin-top: 5px;
}
div.prompt {
display: none;
}
    .CodeMirror{
        font-family: "Consolas", sans-serif;
    }
    
pre, code, kbd, samp {
     font-family: Consolas, monospace;
     font-size: 12px;
}
p {font-size:14px;}
    div.cell{
        max-width:80%;
        margin-left:0px;
        margin-right:auto;
    }
    div.text_cell_render{
        max-width:90%;
        margin-left:0px;
        margin-right:auto;
    }
</style>

<script>
$(document).ready(function(){
  $(".output_wrapper").click(function(){
      $(this).prev('.input_hidden').slideToggle();
  });
})
</script>

{%- endblock header -%}
