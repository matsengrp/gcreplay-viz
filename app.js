// ** Page elements

const url_input = document.getElementById('urlInput');
const my_iframe = document.getElementById('myIframe');
const pdb_inspect_button = document.getElementById('pdbInspectButton');
// const pdb_reset_button = document.getElementById('pdbResetButton');
// const filter_reset_button = document.getElementById('filterResetButton');
const sidebar = document.querySelector('.sidebar');
const sidebar_toggle_button = document.getElementById('sidebarToggleButton');
const alert_message_box = document.getElementById('loadedAlert');
var selector = {
  'pdbid': document.getElementById('pdbidSelector'),
  'metricid': document.getElementById('metricidSelector'),
  'chainid': document.getElementById('chainidSelector'),
};

// ** Data

const github_paths = [
  `https://raw.githubusercontent.com/davidrich27/gcreplay-viz/main`,
]
const github_path = github_paths[0];
const summary_db_path = `${github_path}/data/metadata/summary.json`;
var summary_db = null;

const sidebar_btn_txt = {
  'open': `<<<`,
  'close': `>>>`
};

var prompt = {
  'pdbid': 'Select by PDB',
  'chainid': 'Select by Chain ID',
  'metricid': 'Select by Metric',
};

// Template for dms-viz.github.io query string.
// For query string arguments, see `dms-viz.github.io/v0/js/tool.js:638-696`
var dms_viz_request = {
  'name': null,    /* name displayed at top graph */
  'data': null,    /* url path to json file */
  'pr': 'cartoon', /* focal protein representation */
  'pc': '#b2df8a', /* focal protein chain color */
  'br': 'cartoon', /* peripheral (background) protein representation */
  'bc': '#66ccee', /* peripheral (background) protein chain color */
  'sc': '#ffffff', /* screen color */
  //  'ce': '%5B%22DNSM%22%5D',
};

// ** Functions

class HTMLHelper {
  static raw(text) {
    return `${text}\n`
  }

  static p(text) {
    return `<p>${text}</p>\n`
  }

  static ul_open() {
    return `<ul>\n`
  }

  static ul_close() {
    return `</ul>\n`
  }

  static li(text) {
    return `<li>${text}</li>\n`
  }

  static hr() {
    return `<hr />`
  }

  static ul(text_list) {
    var html_str = HTMLHelper.ul_open();
    for (const i in text_list) {
      const text = text_list[i];
      html_str += HTMLHelper.li(text);
    }
    html_str += HTMLHelper.ul_close();
    return html_str;
  }

  static option_list(opt_values, opt_prompt = null) {
    var select_options = ``
    if (opt_prompt != null) {
      select_options += `\n<option selected>${opt_prompt}: ${opt_values.length}</option>\n`;
    }
    for (const i in opt_values) {
      const value = opt_values[i];
      const select_option = `<option value=\'${value}\'>${value}</option> \n`;
      select_options += select_option;
    }
    return select_options;
  }
}

class Utility {
  static log(text, do_log = true) {
    if (do_log) {
      console.log(`${text}`);
    }
  }

  static load_json(url_path) {
    return new Promise((resolve, reject) => {
      fetch(url_path)
        .then(response => {
          if (!response.ok) {
            throw new Error('Network response was not ok');
          }
          return response.json();
        })
        .then(json_data => {
          resolve(json_data);
        })
        .catch(error => {
          reject(error);
        });
    });

  }

  static async load_summary_db() {
    summary_db = await Utility.load_json(summary_db_path);
    summary_db = new JsonTable(summary_db, 'row');
    return summary_db;
  }

  static iota(length, start = 0, step = 1) {
    const result = [];
    for (let i = 0; i < length; i++) {
      result.push(start + i * step);
    }
    return result;
  }

  static create_unique_sorted_array(arr) {
    const unique_values = new Set(arr);
    const sorted_array = Array.from(unique_values).sort();
    return sorted_array;
  }

  static create_query_string(obj) {
    let query_params = [];
    for (let key in obj) {
      if (obj.hasOwnProperty(key)) {
        query_params.push(encodeURIComponent(key) + '=' + encodeURIComponent(obj[key]));
      }
    }
    return query_params.join('&');
  }

  static truncate(str, delim = '*') {
    return str.split(delim)[0];
  }

  static check_file_exists(file_url, do_log = false) {
    var http = new XMLHttpRequest();
    http.open('HEAD', file_url, false);
    http.send();
    if (http.status === 200) {
      Utility.log(`File exists: ${file_url}`, do_log);
      return true;
    } else {
      Utility.log(`File does NOT exist: ${file_url}`, do_log);
      return false;
    }
  }

  static query_matches(data, query) {
    var is_found = false;
    var split_data = data.split(',')
    for (const i in split_data) {
      const sub_data = split_data[i];
      if (query == sub_data) {
        is_found = true;
        break;
      }
    }
    return is_found;
  }
}

class JsonTable {
  constructor(json_obj, ordered_by = 'col') {
    this.data = json_obj;
    // ordered_by can be ordered first by 'col' or 'row'.
    this.ordered_by = ordered_by;
    this.initialize();
  }

  initialize() {
    this.all_values_by_field = {};
    this.selected_values_by_field = {};
    this.fields = this.get_fields();
    // initialize all fields and selected fields.
    this.all_rowids = this.get_all_rowids();
    this.selected_rowids = this.all_rowids;
    for (const i in this.fields) {
      const field = this.fields[i];
      this.all_values_by_field[field] = this.get_unique_values_by_field(field);
    }
    this.reset_filter();
  }

  get_fields() {
    if (this.ordered_by == 'col') {
      return Object.keys(this.data);
    } /* else if (this.ordered_by == 'row') */
    return Object.keys(this.data[0]);
  }

  get_all_rowids() {
    if (this.ordered_by == 'col') {
      return Object.keys(this.data.pdbid);
    } /* else if (this.orientation == 'row) */
    return Object.keys(this.data);
  }

  get_data(row_id, field_name) {
    if (this.ordered_by == 'col') {
      return this.data[field_name][row_id];
    } /* else if (this.ordered_by == 'row') */
    return this.data[row_id][field_name];
  }

  get_row_data(row_id) {
    var row_data = {};
    for (const i in this.fields) {
      const field_name = this.fields[i];
      row_data[field_name] = this.get_data(row_id, field_name);
    }
    return row_data;
  }

  find_first_row_by_query(field_name, query) {
    for (const i in this.selected_rowids) {
      const row_id = this.selected_rowids[i];
      const data = this.get_data(row_id, field_name);
      if (Utility.query_matches(data, query)) {
        return this.get_row_data(row_id);
      }
    }
    return null;
  }

  reset_filter() {
    this.selected_rowids = this.all_rowids;
    for (const i in this.fields) {
      const field = this.fields[i];
      this.selected_values_by_field[field] = this.all_values_by_field[field];
    }
  }

  update_after_filter() {
    for (const i in this.fields) {
      const field = this.fields[i];
      this.selected_values_by_field[field] = this.get_unique_values_by_field(field);
    }
  }

  apply_filter(field_name, query, update_fields = true) {
    var new_selected_rowids = []
    for (const i in this.selected_rowids) {
      const row_id = this.selected_rowids[i];
      const data = this.get_data(row_id, field_name);
      if (Utility.query_matches(data, query)) {
        new_selected_rowids.push(row_id);
      }
    }
    this.selected_rowids = new_selected_rowids;
    if (update_fields) {
      this.update_after_filter();
    }
    return this.selected_rowids;
  }

  get_values_by_field(field_name) {
    var values = [];
    for (const i in this.selected_rowids) {
      const row_id = this.selected_rowids[i];
      var data = this.get_data(row_id, field_name);
      // data = Utility.truncate(data, '*');
      values.push(data);
    }
    return values;
  }

  get_unique_values_by_field(field_name) {
    var vals = this.get_values_by_field(field_name);
    return Utility.create_unique_sorted_array(vals);
  }
}

class Event {
  // Window slider and alerts

  static sidebar_toggle() {
    sidebar.classList.toggle('sidebar-collapse');
    sidebar_toggle_button.innerText = (sidebar_toggle_button.innerText == sidebar_btn_txt['open']) ? sidebar_btn_txt['close'] : sidebar_btn_txt['open'];
  }

  static alert_set_text(text, append = false, alert_color = null) {
    if (append) {
      alert_message_box.innerHTML += HTMLHelper.hr();
    } else {
      alert_message_box.innerHTML = '';
    }
    alert_message_box.innerHTML += HTMLHelper.p(text);

    if (alert_color) {
      Event.alert_color(alert_color);
    }
  }

  static alert_set_color(alert_color) {
    if (alert_color != null) {
      alert_message_box.classList.remove('alert-danger');
      alert_message_box.classList.remove('alert-warning');
      alert_message_box.classList.remove('alert-success');
      alert_message_box.classList.add(alert_color);
    }
  }

  // Submitting requests

  static submit_dms_viz_request(request = dms_viz_request) {
    var query_string = Utility.create_query_string(request);
    var url = `https://dms-viz.github.io/v0/?${query_string}`;
    console.log(`Loading url into iframe: ${url}`);
    my_iframe.src = url;
  }

  static load_pdb() {
    // get field from database
    const filter = {};
    for (const key in selector) {
      const value = selector[key].value;
      if (value !== "") {
        filter[key] = value;
      }
    }
    const matches = summary_db.data.filter(item =>
      Object.entries(filter).every(([key, val]) => item[key] == val)
    );
    const match = matches[0];
    console.log(`matches found: ${match.length}`)
    console.log(matches)

    // get copy of request and fill fields
    const my_dms_viz_request = JSON.parse(JSON.stringify(dms_viz_request));
    my_dms_viz_request.name = match['description'];
    const file_name = match['dmsviz_filepath'];
    my_dms_viz_request.data = `${github_path}/data/dmsviz-jsons/${file_name}`;

    // log request
    var alert_text = `Loading pdb...\n`
    alert_text += `${HTMLHelper.ul_open()}`
    alert_text += `<li> PDB: ${match['pdbid']} </li>\n`
    alert_text += `<li> Chain ID: ${match['chainid']} </li>\n`
    alert_text += `<li> Metric: ${match['metric_full_name']} </li>\n`
    alert_text += `${HTMLHelper.ul_close()}`
    Event.alert_set_text(alert_text)

    console.log(my_dms_viz_request);
    Event.submit_dms_viz_request(my_dms_viz_request);
  }

  static load_pdb_from_query_string() {
    const url = new URL(window.location.href);
    const params = new URLSearchParams(url.search.toLowerCase());
    console.log(params);
    const query_pdbid = params.get('pdbid');
    const query_metric = params.get('metric')
    var is_pdb_loaded = false;
    if (query_pdbid) {
      is_pdb_loaded = Event.load_pdb(query_pdbid);
      if (!is_pdb_loaded) {
        alert_message_box.innerText = `PDBID requested from querystring does not exist: ${query_pdbid}.`;
        alert_message_box.classList.add('alert-danger');
      }
    }
    return is_pdb_loaded;
  }

  // Menus

  static get_first_unique_rows_by_column(data, field) {
    const seen = new Map();
    data.forEach(row => {
      const key = row[field];
      if (!seen.has(key)) {
        seen.set(key, row);
      }
    });
    return Array.from(seen.values());
  }

  static populate_dropdown_from_data(menu_elem, data, value_field, text_field) {
    const prompt_text = prompt[value_field]
    menu_elem.innerHTML = `<option value="">-- ${prompt_text} --</option>`;

    // const unique_values = [...new Set(data.map(item => item[value_field]))].sort();
    const unique_rows = Event.get_first_unique_rows_by_column(data, value_field);

    unique_rows.forEach(item => {
      const option = document.createElement('option');
      option.value = item[value_field];
      option.text = item[text_field];
      menu_elem.appendChild(option);
    });
  }
}

// ** Event Listeners

document.addEventListener('DOMContentLoaded', async function () {
  // Initialize page elements
  sidebar_toggle_button.innerText = sidebar_btn_txt['open']

  // Initialize data
  summary_db = await Utility.load_summary_db();
  // summary_db = new JsonTable(summary_db, 'row');

  // Populate data
  Event.populate_dropdown_from_data(selector['pdbid'], summary_db.data, 'pdbid', 'pdbid');
  Event.populate_dropdown_from_data(selector['chainid'], summary_db.data, 'chainid', 'chainid');
  Event.populate_dropdown_from_data(selector['metricid'], summary_db.data, 'metricid', 'metric_full_name');

  // Event buttons
  sidebar_toggle_button.addEventListener('click', Event.sidebar_toggle);
  pdb_inspect_button.addEventListener('click', Event.load_pdb);
});
