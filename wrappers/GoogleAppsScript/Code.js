function PropsSI(output,property_1,value_1,property_2,value_2,fluid, opt_url) {

  /* This function fetches the state property for a fluid from the web api

  Args:
  output: str, name of the output property
  property_1: str, name of the first input property
  value_1: float, value of the first input property
  property_2: str, name of the second input property
  value_2: float, value of the second input property
  fluid: str, name of the fluid
  opt_url: str optional, url of the web api

  Returns:
  float, value of the state property
  
  */

  if (opt_url==null) {
    var baseUrl = "https://qnikil7.pythonanywhere.com"
  } else {
    var baseUrl = opt_url
  }
  
  functionName = "PropsSI"

  var requestString = baseUrl + '/' + functionName + '/' + output + '/' + property_1 + '/' + value_1 + '/' + property_2 + '/' + value_2 + '/' + fluid ;

  var response = UrlFetchApp.fetch(requestString, {muteHttpExceptions: true});

  var responseData = JSON.parse(response.getContentText());

  return responseData;
  
}

function HAPropsSI(output,property_1,value_1,property_2,value_2, property_3,value_3, opt_url) {

  /* This function fetches the humid air state property from the web api

  Args:
  output: str, name of the output property
  property_1: str, name of the first input property
  value_1: float, value of the first input property
  property_2: str, name of the second input property
  value_2: float, value of the second input property
  property_3: str, name of the third input property
  value_3: float, value of the third input property
  opt_url: str optional, url of the web api

  Returns:
  float, value of the state property

  */

  if (opt_url==null) {
    var baseUrl = "https://qnikil7.pythonanywhere.com"
  } else {
    var baseUrl = opt_url
  }
  
  functionName = "HAPropsSI"

  var requestString = baseUrl + '/' + functionName + '/' + output + '/' + property_1 + '/' + value_1 + '/' + property_2 + '/' + value_2 + '/' + property_3 + '/' + value_3 ;

  var response = UrlFetchApp.fetch(requestString, {muteHttpExceptions: true});

  var responseData = JSON.parse(response.getContentText());

  return responseData;
  
}
