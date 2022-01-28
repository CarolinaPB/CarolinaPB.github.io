---
layout: post
title: Add Snakemake snippets to VS Code
tags: [snakemake, Virtual Studio Code, snippets, bioinformatics]
---

Writting a Snakemake pipeline can be quite repetitive, you add rule after rule and they all follow the same format. To make it easier and faster to write a Snakemake pipeline I created [this snippets file](https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/snakemake_rule.code-snippets) in Virtual Studio Code (it may also be possible in other editors).    
To create a snippets file in VS Code you go to `file` -> `Preferences` -> `User Snippets`, then you choose your language, in this case `snakemake`. The snippets that you create will be available when you have a snakemake file open.  
Once you have created the snippets, to use them you just need to start typing the prefix and hit `enter`.

These are the snippets I've been using:
## 1. Add rule body
```
"Add rule": {
  "scope": "snakemake",
  "prefix": "rule",
  "body": [
    "rule ${1:NAME}$0:",
    "\tinput:",
    "\t\t'${2:INPUT}'",
    "\toutput:",
    "\t\t'${3:OUTPUT}'",
    "\tmessage:",
    "\t\t'Rule {rule} processing'",
    "\tshell:",
    "\t\t'${4}'"
  ],
  "description": "Add rule body"
},
```
With this snippet, you start typing "rule" and VSC will suggest the snippet. The only thing you need to do is press `enter`.  
A new rule will be created:
```
rule NAME:
    input:
        'INPUT'
    output:
        'OUTPUT'
    message:
        'Rule {rule} processing'
    shell:
        ''
```
The rule has `input`, `output`, `message`, and `shell` fields. It also has placeholders: `NAME`, `INPUT`, `OUTPUT`, and one in the shell `''`.   
You can start editing your rule right away. The cursor will be at the end of `NAME` and you can immediately start typing to change the rule's name. You can use `TAB` to move to the next placeholder. 

## 2. Add message
```
"Add message to rule":{
  "scope": "snakemake",
  "prefix": "message",
  "body": [
    "message:",
    "\t'Rule {rule} processing'"
  ],
  "description": "Add message component to rule"
},
```
This snippet can be used by starting to write `message` and pressing `enter`.  
Use this snippet to add a `message` field to your existing rule:
```
message:
    'Rule {rule} processing'
```

## 3. Add group

```
"Add group to rule": {
  "scope": "snakemake",
  "prefix": "group",
  "body": [
    "group:",
    "\t'group'"
  ],
  "description": "Add group component to rule"
},
```
This snippet can be used by starting to write `group` and pressing `enter`.  
Use this snippet to add a `group` field to your existing rule:

```
group:
    'group'
```

### Complete file
```
{
	// Place your global snippets here. Each snippet is defined under a snippet name and has a scope, prefix, body and 
	// description. Add comma separated ids of the languages where the snippet is applicable in the scope field. If scope 
	// is left empty or omitted, the snippet gets applied to all languages. The prefix is what is 
	// used to trigger the snippet and the body will be expanded and inserted. Possible variables are: 
	// $1, $2 for tab stops, $0 for the final cursor position, and ${1:label}, ${2:another} for placeholders. 
	// Placeholders with the same ids are connected.
	// Example:
	// "Print to console": {
	// 	"scope": "javascript,typescript",
	// 	"prefix": "log",
	// 	"body": [
	// 		"console.log('$1');",
	// 		"$2"
	// 	],
	// 	"description": "Log output to console"
	// }
	
	"Add rule": {
		"scope": "snakemake",
		"prefix": "rule",
		"body": [
			"rule ${1:NAME}$0:",
			"\tinput:",
			"\t\t'${2:INPUT}'",
			"\toutput:",
			"\t\t'${3:OUTPUT}'",
			"\tmessage:",
			"\t\t'Rule {rule} processing'",
			"\tshell:",
			"\t\t'${4}'"
		],
		"description": "Add rule body"
	},
	
	"Add message to rule":{
		"scope": "snakemake",
		"prefix": "message",
		"body": [
			"message:",
			"\t'Rule {rule} processing'"
		],
		"description": "Add message component to rule"
	},

	"Add group to rule": {
		"scope": "snakemake",
		"prefix": "group",
		"body": [
			"group:",
			"\t'group'"
		],
		"description": "Add group component to rule"
	},

}
```

To read more about snippets in VS Code: https://code.visualstudio.com/docs/editor/userdefinedsnippets
